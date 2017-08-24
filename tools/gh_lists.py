#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
"""gh_lists.py [options] MILESTONE

Extract pull request and issue lists for a given milestone from Github.

Also insert text fragments to include in the release notes.


Release note fragments
----------------------

Release note fragments can be included in Github pull requests in
format:

[release-notes type="feature" section="scipy.cluster"]
Text to insert, in Github markdown format. For major new features,
include a 3rd level section header (`### some title`) on top.
[/release-notes]

The "type" can be one of "feature" (new improvement), "change"
(backward incompatible change), "deprecation" (deprecated message),
and "other" (any other change).

There may be multiple such sections in each pull request.

Github markdown will be converted to restructuredtext using pandoc,
which needs to be installed.

"""
from __future__ import print_function, division, absolute_import

import os
import re
import sys
import json
import time
import tempfile
import datetime
import textwrap
import collections
import subprocess
import argparse

from urllib.request import urlopen, Request, HTTPError
from tempita import Template

TEMPLATE_FILE = os.path.join(os.path.dirname(__file__), 'gh_lists.rst')

Issue = collections.namedtuple('Issue', ('id', 'title', 'url'))
ReleaseNoteFragment = collections.namedtuple('ReleaseNoteFragment',
                                             ('issue_id', 'type', 'section', 'text',
                                              'has_title'))


def main():
    p = argparse.ArgumentParser(usage=__doc__.lstrip())
    p.add_argument('--project', default='scipy/scipy',
                   help="Github project to fetch milestone information from")
    p.add_argument('--auth', action='store_true',
                   help="authenticate to Github (increases rate limits)")
    p.add_argument('milestone')
    args = p.parse_args()

    cache_file = os.path.join(os.path.dirname(__file__), 'gh_cache.json')
    getter = GithubGet(cache_file, auth=args.auth)
    try:
        milestones = get_milestones(getter, args.project)
        if args.milestone not in milestones:
            msg = "Milestone {0} not available. Available milestones: {1}"
            msg = msg.format(args.milestone, u", ".join(sorted(milestones)))
            p.error(msg)
        issues, fragments = get_issues(getter, args.project, args.milestone)
        issues.sort()
        fragments.sort(key=lambda f: (not f.has_title, f))
    finally:
        getter.save()

    prs = [x for x in issues if u'/pull/' in x.url]
    issues = [x for x in issues if x not in prs]

    output = format_template(args.milestone, prs, issues, fragments)

    print(output)
    return 0


def format_template(milestone, prs, issues, fragments):
    template = Template.from_filename(TEMPLATE_FILE, encoding='utf-8')

    def warn(message):
        print("WARNING:", message, file=sys.stderr, flush=True)

    def format_title(title, underline="-"):
        title = title.strip()
        return title.strip() + "\n" + underline*len(title) + "\n"

    def format_list(items):
        parts = []
        for issue in items:
            msg = u"- `#{0} <{1}>`__: {2}"
            title = re.sub(u"\s+", u" ", issue.title.strip())
            if len(title) > 60:
                remainder = re.sub(u"\s.*$", u"...", title[60:])
                if len(remainder) > 20:
                    remainder = title[:80] + u"..."
                else:
                    title = title[:60] + remainder
            parts.append(msg.format(issue.id, issue.url, title))
        return "\n".join(parts)

    types = sorted(set(fragment.type for fragment in fragments))

    sections = {}
    for typ in types:
        secs = sorted(set(fragment.section for fragment in fragments
                          if fragment.type == typ))
        sections[typ] = {}
        for sec in secs:
            sections[typ][sec] = sorted(fragment for fragment in fragments
                                        if fragment.type == typ and fragment.section == sec)

    context = dict(milestone=milestone,
                   issues=issues,
                   prs=prs,
                   all_fragments=fragments,
                   sections=sections,
                   types=types,
                   warn=warn,
                   format_title=format_title,
                   format_list=format_list)

    return template.substitute(context).strip()


def get_milestones(getter, project):
    url = "https://api.github.com/repos/{project}/milestones".format(project=project)
    data, info = getter.get(url)

    milestones = {}
    for ms in data:
        milestones[ms[u'title']] = ms[u'number']
    return milestones


def get_issues(getter, project, milestone):
    milestones = get_milestones(getter, project)
    mid = milestones[milestone]

    url = "https://api.github.com/repos/{project}/issues?milestone={mid}&state=closed&sort=created&direction=asc"
    url = url.format(project=project, mid=mid)

    data = getter.get_multipage(url)

    issues = []
    fragments = []

    for issue_data in data:
        issues.append(Issue(issue_data[u'number'],
                            issue_data[u'title'],
                            issue_data[u'html_url']))
        fragments.extend(parse_fragments(issue_data[u'number'],
                                         issue_data[u'body']))

    return issues, fragments


def parse_fragments(number, text):
    if not text:
        return []

    items = re.findall(r'\[release-notes([^\]]*)\](.*?)\[/release-notes\]',
                       text, flags=re.I | re.S)

    fragments = []
    for flags, text_md in items:
        m = re.search('section\\s*=\\s*["\'](.*?)["\']', flags)
        if m:
            section = m.group(1).strip().lower()
        else:
            section = None

        m = re.search('type\\s*=\\s*["\'](.*?)["\']', flags)
        if m:
            typ = m.group(1).strip().lower()
        else:
            typ = None

        text_md = text_md.strip()

        if text_md.startswith('#'):
            # Ensure correct title level
            text_md = re.sub('^#+', '###', text_md, flags=re.S)
            has_title = True
        else:
            has_title = False

        p = subprocess.Popen(['pandoc', '-f', 'markdown', '-t', 'rst'], stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)
        out, err = p.communicate(text_md.encode('utf-8'))
        if p.returncode != 0:
            print("WARNING: gh-{}: pandoc failed to convert markdown to rst: {}\n{}".format(number, err, out),
                  file=sys.stderr, flush=True)
            print("Source:\n{}".format(text_md), file=sys.stderr, flush=True)
            text_rst = text_md
        else:
            text_rst = out.decode('utf-8')

        fragments.append(ReleaseNoteFragment(number, typ, section, text_rst,
                                             has_title))

    return fragments


class GithubGet(object):
    def __init__(self, filename, auth=False):
        self.headers = {'User-Agent': 'scipy/gh_lists.py',
                        'Accept': 'application/vnd.github.v3+json'}

        if auth:
            self.authenticate()

        req = self.urlopen('https://api.github.com/rate_limit')
        try:
            if req.getcode() != 200:
                raise RuntimeError()
            info = json.loads(req.read().decode('utf-8'))
        finally:
            req.close()

        self.ratelimit_remaining = int(info['rate']['remaining'])
        self.ratelimit_reset = float(info['rate']['reset'])

        self.filename = filename
        if os.path.isfile(filename):
            print("[gh_lists] using {0} as cache (remove it if you want fresh data)".format(filename),
                  file=sys.stderr)
            with open(filename, 'r', encoding='utf-8') as f:
                self.cache = json.load(f)
        else:
            self.cache = {}

    def authenticate(self):
        print("Input a Github API access token.\n"
              "Personal tokens can be created at https://github.com/settings/tokens\n"
              "This script does not require any permissions (so don't give it any).",
              file=sys.stderr, flush=True)
        print("Access token: ", file=sys.stderr, end='', flush=True)
        token = input()
        self.headers['Authorization'] = 'token {0}'.format(token.strip())

    def urlopen(self, url, auth=None):
        assert url.startswith('https://')
        req = Request(url, headers=self.headers)
        return urlopen(req, timeout=60)

    def get_multipage(self, url):
        data = []
        while url:
            page_data, info = self.get(url)
            data += page_data
            url = info['Next']
        return data

    def _get(self, url):
        while True:
            # Wait until rate limit
            while self.ratelimit_remaining == 0 and self.ratelimit_reset > time.time():
                s = self.ratelimit_reset + 5 - time.time()
                if s <= 0:
                    break
                print("[gh_lists] rate limit exceeded: waiting until {0} ({1} s remaining)".format(
                         datetime.datetime.fromtimestamp(self.ratelimit_reset).strftime('%Y-%m-%d %H:%M:%S'),
                         int(s)),
                      file=sys.stderr, flush=True)
                time.sleep(min(5*60, s))

            # Get page
            print("[gh_lists] get:", url, file=sys.stderr, flush=True)
            try:
                req = self.urlopen(url)
                try:
                    code = req.getcode()
                    info = dict(req.info())
                    data = json.loads(req.read().decode('utf-8'))
                finally:
                    req.close()
            except HTTPError as err:
                code = err.getcode()
                info = err.info()
                data = None

            if code not in (200, 403):
                raise RuntimeError()

            # Parse reply
            info['Next'] = None
            if 'Link' in info:
                m = re.search('<([^<>]*)>; rel="next"', info['Link'])
                if m:
                    info['Next'] = m.group(1)

            # Update rate limit info
            if 'X-RateLimit-Remaining' in info:
                self.ratelimit_remaining = int(info['X-RateLimit-Remaining'])
            if 'X-RateLimit-Reset' in info:
                self.ratelimit_reset = float(info['X-RateLimit-Reset'])

            # Deal with rate limit exceeded
            if code != 200 or data is None:
                if self.ratelimit_remaining == 0:
                    continue
                else:
                    raise RuntimeError()

            # Done.
            return data, info

    def get(self, url):
        if url not in self.cache:
            self.cache[url] = self._get(url)
        else:
            print("[gh_lists] get (cached):", url, file=sys.stderr)
        return self.cache[url]

    def save(self):
        f = tempfile.NamedTemporaryFile(mode='w',
                                        encoding='utf-8',
                                        dir=os.path.dirname(self.filename),
                                        prefix=os.path.basename(self.filename) + '-',
                                        delete=False)
        try:
            json.dump(self.cache, f)
            f.close()
            os.rename(f.name, self.filename)
        except:
            if not f.closed:
                f.close()
            os.unlink(f.name)
            raise


def test_parse_fragments():
    body = textwrap.dedent("""
    Some text
    Some more text

    [release-notes type="feature" section="scipy.sparse"]
    # Some big feature

    We like major features, and want more of them.
    See links at [scipy.org](https://www.scipy.org).
    [/release-notes]
    """)
    expected = textwrap.dedent("""\
    Some big feature
    ~~~~~~~~~~~~~~~~

    We like major features, and want more of them. See links at
    `scipy.org <https://www.scipy.org>`__.
    """)

    fragments = parse_fragments(123, body)
    assert fragments == [ReleaseNoteFragment(section="scipy.sparse",
                                             type="feature",
                                             text=expected,
                                             issue_id=123,
                                             has_title=True)]

    result = format_template("1.0.0", [], [], fragments)
    assert 'Some big feature' in result


if __name__ == "__main__":
    sys.exit(main())

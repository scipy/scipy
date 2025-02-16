#!/usr/bin/env python3
"""
gh_lists.py MILESTONE

Functions for Github API requests.
"""
import os
import re
import sys
import json
import collections
import argparse
import datetime
import time

from urllib.request import urlopen, Request, HTTPError


Issue = collections.namedtuple('Issue', ('id', 'title', 'url'))


def main():
    p = argparse.ArgumentParser(usage=__doc__.lstrip())
    p.add_argument('--project', default='scipy/scipy')
    p.add_argument('milestone')
    args = p.parse_args()

    getter = CachedGet('gh_cache.json', GithubGet())
    try:
        milestones = get_milestones(getter, args.project)
        if args.milestone not in milestones:
            msg = "Milestone {0} not available. Available milestones: {1}"
            msg = msg.format(args.milestone, ", ".join(sorted(milestones)))
            p.error(msg)
        issues = get_issues(getter, args.project, args.milestone)
        issues.sort()
    finally:
        getter.save()

    prs = [x for x in issues if '/pull/' in x.url]
    issues = [x for x in issues if x not in prs]

    def print_list(title, items):
        print()
        print(title)
        print("-"*len(title))
        print()

        def backtick_repl(matchobj):
            """repl to add an escaped space following a code block if needed"""
            if matchobj.group(2) != ' ':
                post = r'\ ' + matchobj.group(2)
            else:
                post = matchobj.group(2)
            return '``' + matchobj.group(1) + '``' + post
        
        def asterisk_repl(matchobj):
            """repl to un-escape asterisks in code blocks"""
            code = matchobj.group(1).replace("\\*", "*")
            return '``' + code + '``'

        for issue in items:
            msg = "* `#{0} <{1}>`__: {2}"
            # sanitize whitespace
            title = re.sub("\\s+", " ", issue.title.strip())

            # substitute any single backtick not adjacent to a backtick
            # for a double backtick
            title = re.sub("(?P<pre>(?:^|(?<=[^`])))`(?P<post>(?=[^`]|$))",
                           r"\g<pre>``\g<post>",
                           title)
            # add an escaped space if code block is not followed by a space
            title = re.sub("``(.*?)``(.)", backtick_repl, title)

            # sanitize asterisks
            title = title.replace('*', '\\*')
            # except those in code blocks
            title = re.sub("``(.*?)``", asterisk_repl, title)

            if len(title) > 60:
                remainder = re.sub("\\s.*$", "...", title[60:])
                if len(remainder) > 20:
                    # just use the first 80 characters, with ellipses.
                    # note: this was previously bugged,
                    # assigning to `remainder` rather than `title`
                    title = title[:80] + "..."
                else:
                    # use the first 60 characters and the next word
                    title = title[:60] + remainder

                if title.count('`') % 4 != 0:
                    # ellipses have cut in the middle of a code block,
                    # so finish the code block before the ellipses
                    title = title[:-3] + '``...'

            msg = msg.format(issue.id, issue.url, title)
            print(msg)
        print()

    msg = f"Issues closed for {args.milestone}"
    print_list(msg, issues)

    msg = f"Pull requests for {args.milestone}"
    print_list(msg, prs)

    return 0


def get_milestones(getter, project):
    url = f"https://api.github.com/repos/{project}/milestones"
    data = getter.get(url)

    milestones = {}
    for ms in data:
        milestones[ms['title']] = ms['number']
    return milestones


def get_issues(getter, project, milestone):
    milestones = get_milestones(getter, project)
    mid = milestones[milestone]

    url = "https://api.github.com/repos/{project}/issues?milestone={mid}&state=closed&sort=created&direction=asc"
    url = url.format(project=project, mid=mid)

    data = getter.get(url)

    issues = []
    # data contains both PR and issue data
    for issue_data in data:
        # don't include PRs that were closed instead
        # of merged
        if "pull" in issue_data['html_url']:
            merge_status = issue_data['pull_request']['merged_at']
            if merge_status is None:
                continue
        issues.append(Issue(issue_data['number'],
                            issue_data['title'],
                            issue_data['html_url']))
    return issues


class CachedGet:
    def __init__(self, filename, getter):
        self._getter = getter

        self.filename = filename
        if os.path.isfile(filename):
            print(f"[gh_lists] using {filename} as cache "
                  f"(remove it if you want fresh data)",
                  file=sys.stderr)
            with open(filename, encoding='utf-8') as f:
                self.cache = json.load(f)
        else:
            self.cache = {}

    def get(self, url):
        if url not in self.cache:
            data = self._getter.get_multipage(url)
            self.cache[url] = data
            return data
        else:
            print("[gh_lists] (cached):", url, file=sys.stderr, flush=True)
            return self.cache[url]

    def save(self):
        tmp = self.filename + ".new"
        with open(tmp, 'w', encoding='utf-8') as f:
            json.dump(self.cache, f)
        os.rename(tmp, self.filename)


class GithubGet:
    def __init__(self, auth=False):
        self.headers = {'User-Agent': 'gh_lists.py',
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

    def authenticate(self):
        print("Input a Github API access token.\n"
              "Personal tokens can be created at https://github.com/settings/tokens\n"
              "This script does not require any permissions (so don't give it any).",
              file=sys.stderr, flush=True)
        print("Access token: ", file=sys.stderr, end='', flush=True)
        token = input()
        self.headers['Authorization'] = f'token {token.strip()}'

    def urlopen(self, url, auth=None):
        assert url.startswith('https://')
        req = Request(url, headers=self.headers)
        return urlopen(req, timeout=60)

    def get_multipage(self, url):
        data = []
        while url:
            page_data, info, next_url = self.get(url)
            data += page_data
            url = next_url
        return data

    def get(self, url):
        while True:
            # Wait until rate limit
            while self.ratelimit_remaining == 0 and self.ratelimit_reset > time.time():
                s = self.ratelimit_reset + 5 - time.time()
                if s <= 0:
                    break
                print(
                    "[gh_lists] rate limit exceeded: waiting until {} ({} s remaining)"
                    .format(datetime.datetime.fromtimestamp(self.ratelimit_reset)
                            .strftime('%Y-%m-%d %H:%M:%S'),
                            int(s)),
                    file=sys.stderr, flush=True
                )
                time.sleep(min(5*60, s))

            # Get page
            print("[gh_lists] get:", url, file=sys.stderr, flush=True)
            try:
                req = self.urlopen(url)
                try:
                    code = req.getcode()
                    info = req.info()
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
            next_url = None
            if 'Link' in info:
                m = re.search('<([^<>]*)>; rel="next"', info['Link'])
                if m:
                    next_url = m.group(1)

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
            return data, info, next_url


if __name__ == "__main__":
    sys.exit(main())

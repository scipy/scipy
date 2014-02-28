#!/usr/bin/env python
# -*- encoding:utf-8 -*-
"""
gh_lists.py MILESTONE

Functions for Github API requests.
"""
from __future__ import print_function, division, absolute_import

import os
import re
import sys
import json
import collections
import argparse

from urllib2 import urlopen


Issue = collections.namedtuple('Issue', ('id', 'title', 'url'))


def main():
    p = argparse.ArgumentParser(usage=__doc__.lstrip())
    p.add_argument('--project', default='scipy/scipy')
    p.add_argument('milestone')
    args = p.parse_args()

    getter = CachedGet('gh_cache.json')
    try:
        milestones = get_milestones(getter, args.project)
        if args.milestone not in milestones:
            msg = "Milestone {0} not available. Available milestones: {1}"
            msg = msg.format(args.milestone, u", ".join(sorted(milestones)))
            p.error(msg)
        issues = get_issues(getter, args.project, args.milestone)
        issues.sort()
    finally:
        getter.save()

    prs = [x for x in issues if u'/pull/' in x.url]
    issues = [x for x in issues if x not in prs]

    def print_list(title, items):
        print()
        print(title)
        print("-"*len(title))
        print()

        for issue in items:
            msg = u"- `#{0} <{1}>`__: {2}"
            title = re.sub(u"\s+", u" ", issue.title.strip())
            if len(title) > 60:
                remainder = re.sub(u"\s.*$", u"...", title[60:])
                if len(remainder) > 20:
                    remainder = title[:80] + u"..."
                else:
                    title = title[:60] + remainder
            msg = msg.format(issue.id, issue.url, title)
            print(msg)
        print()

    msg = u"Issues closed for {0}".format(args.milestone)
    print_list(msg, issues)
    
    msg = u"Pull requests for {0}".format(args.milestone)
    print_list(msg, prs)

    return 0


def get_milestones(getter, project):
    url = "https://api.github.com/repos/{project}/milestones".format(project=project)
    raw_data, info = getter.get(url)
    data = json.loads(raw_data)

    milestones = {}
    for ms in data:
        milestones[ms[u'title']] = ms[u'number']
    return milestones


def get_issues(getter, project, milestone):
    milestones = get_milestones(getter, project)
    mid = milestones[milestone]

    url = "https://api.github.com/repos/{project}/issues?milestone={mid}&state=closed&sort=created&direction=asc"
    url = url.format(project=project, mid=mid)

    raw_datas = []
    while True:
        raw_data, info = getter.get(url)
        raw_datas.append(raw_data)
        if 'link' not in info:
            break
        m = re.search('<(.*?)>; rel="next"', info['link'])
        if m:
            url = m.group(1)
            continue
        break

    issues = []

    for raw_data in raw_datas:
        data = json.loads(raw_data)
        for issue_data in data:
            issues.append(Issue(issue_data[u'number'],
                                issue_data[u'title'],
                                issue_data[u'html_url']))
    return issues


class CachedGet(object):
    def __init__(self, filename):
        self.filename = filename
        if os.path.isfile(filename):
            print("[gh_lists] using {0} as cache (remove it if you want fresh data)".format(filename),
                  file=sys.stderr)
            with open(filename, 'rb') as f:
                self.cache = json.load(f)
        else:
            self.cache = {}

    def get(self, url):
        url = unicode(url)
        if url not in self.cache:
            print("[gh_lists] get:", url, file=sys.stderr)
            req = urlopen(url)
            if req.getcode() != 200:
                raise RuntimeError()
            data = req.read()
            info = dict(req.info())
            self.cache[url] = (data, info)
            req.close()
        else:
            print("[gh_lists] get (cached):", url, file=sys.stderr)
        return self.cache[url]

    def save(self):
        tmp = self.filename + ".new"
        with open(tmp, 'wb') as f:
            json.dump(self.cache, f)
        os.rename(tmp, self.filename)


if __name__ == "__main__":
    sys.exit(main())

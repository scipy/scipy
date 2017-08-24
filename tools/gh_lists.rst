{{# gh_lists.rst -- Tempita template for gh_lists.py output

The following items are predefined in the context:

milestone = "X.Y.Z"
issues = list of issues: [namedtuple('Issue', ('id', 'title', 'url')), ...]
prs = list of pull requests: [namedtuple('Issue', ('id', 'title', 'url')), ...]
all_fragments = release note fragments: namedtuple('ReleaseNoteFragment',
                                             ('issue_id', 'type', 'section', 'text',
                                              'has_title'))
sections = {type: {section: [fragment1, fragment2, ...], ...}, ...}
types = list of types of fragments present

def warn(msg): emit warning message
def format_title(title, underline="-"): return rst underlined title
def format_list(issue_list): return rst formatted list of issues
}}

{{py:

TYPES = ["feature", "deprecation", "change", "other"]

new_fragments = []
for fragment in all_fragments:
    if fragment.type not in TYPES and fragment.type != "other":
        warn("gh-{}: unknown fragment type {}".format(fragment.issue_id, fragment.type))
}}

==========================
SciPy {{milestone}} Release Notes
==========================

{{if sections.get("feature")}}
New features
============


{{for section, fragments in sorted(sections["feature"].items())}}
{{format_title("`{}` improvements".format(section))}}
{{for fragment in fragments}}
{{fragment.text}}
{{endfor}}
{{endfor}}
{{endif}}

{{if sections.get("deprecated")}}
Deprecated features
===================


{{for section, fragments in sorted(sections["deprecation"].items())}}
{{for fragment in fragments}}
{{fragment.text}}
{{endfor}}
{{endfor}}
{{endif}}

{{if sections.get("change")}}
Backwards incompatible changes
==============================


{{for section, fragments in sorted(sections["change"].items())}}
{{for fragment in fragments}}
{{fragment.text}}
{{endfor}}
{{endfor}}
{{endif}}

Other changes
=============


{{for types, secs in sorted(sections.items())}}
{{for section, fragments in sorted(secs.items())}}
{{for fragment in fragments}}
{{if fragment.type not in ["feature", "deprecation", "change"]}}
{{fragment.text}}
{{endif}}
{{endfor}}
{{endfor}}
{{endfor}}


{{format_title("Issues closed for {}".format(milestone), "=")}}
{{format_list(issues)}}

{{format_title("Pull requests for {}".format(milestone), "=")}}
{{format_list(prs)}}

.. _triaging:

============================
Triaging and curating issues
============================

SciPy has many hundreds of open issues. Closing invalid ones and correctly
labelling valid ones (ideally with some first thoughts in a comment) allows
prioritizing maintenance work and finding related issues easily when working on
an existing function or subpackage.

While anyone can comment and give more information on open issues, extra
permissions are needed if you want to apply labels to issues in the SciPy
repository. While there is no formal process to receive triage rights, the
expectation is that someone should be active as a contributor before joining
the team.

Roles and permissions
=====================

SciPy uses two levels of permissions: triage and core team members. **Triage
members** can label and close issues and pull requests, while **maintainers**
can label and close issues and pull request, and can also merge pull requests.

`GitHub publishes the full list of permissions for the platform.
<https://docs.github.com/en/organizations/managing-access-to-your-organizations-repositories/repository-roles-for-an-organization>`__

Improving issues
================

Issue descriptions can be incomplete, inaccurate or outdated. No special
permissions are needed to work on improving them - this can be useful and help
reduce the workload for maintainers and other contributors. The following
actions are typically useful:

- documenting issues that are missing elements to reproduce the problem such as
  code samples
- suggesting to reformulate the title and description to make them more explicit
  about the problem to be solved
- linking to related issues or discussions while briefly describing how they are
  related, for instance “See also #xyz for a similar attempt at this” provides
  context and helps the discussion.

Keep in mind that every comment on an issue or pull request creates a
notification for a group of people. Be mindful and make use of the edit comment
button when necessary.

Fruitful discussions
====================

Online discussions may be harder than it seems at first glance, in particular
given that a person new to open-source may have a very different understanding
of the process than a seasoned maintainer. 

Overall, it is useful to stay positive and assume good will.
`This article <http://gael-varoquaux.info/programming/technical-discussions-are-hard-a-few-tips.html>`__
explores how to lead online discussions in the context of open source. It's also
important to remember that all interactions are expected to follow the
:ref:`SciPy Code of Conduct <scipy-coc>`.

Issue labels (requires triage rights)
=====================================

When an issue or pull request is created, SciPy may automatically assign one or
more labels depending on the title or section of the code involved. For example,
all issues created with title including the ``BUG:`` prefix will automatically
receive a ``defect`` label.

In some cases, it may be useful to also include other labels manually. Any
person with triage rights can add or remove label as appropriate. Check the
`full description of the current labels <https://github.com/scipy/scipy/labels>`__
for more information.

Other references
================

* `scikit-learn's documentation on Bug triaging and issue curation <https://scikit-learn.org/dev/developers/bug_triaging.html>`__
* `pandas documentation on Issue triage <https://pandas.pydata.org/docs/development/maintaining.html>`__

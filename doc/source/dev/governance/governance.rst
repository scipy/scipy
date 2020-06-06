.. _governance:

========================
SciPy Project Governance
========================

The purpose of this document is to formalize the governance process
used by the SciPy project in both ordinary and extraordinary
situations, and to clarify how decisions are made and how the various
elements of our community interact, including the relationship between
open source collaborative development and work that may be funded by
for-profit or non-profit entities.


The Project
===========

The SciPy Project (The Project) is an open source software project.
The goal of The Project is to develop open source software for scientific
computing in Python, and, in particular, the ``scipy`` package. The Software
developed by The Project is released under the BSD (or similar) open source
license, developed openly and hosted on public GitHub repositories under
the ``scipy`` GitHub organization.

The Project is developed by a team of distributed developers, called
Contributors. Contributors are individuals who have contributed code,
documentation, designs, or other work to the Project. Anyone can be a
Contributor. Contributors can be affiliated with any legal entity or
none. Contributors participate in the project by submitting, reviewing,
and discussing GitHub Pull Requests and Issues and participating in open
and public Project discussions on GitHub, mailing lists, and other
channels. The foundation of Project participation is openness and
transparency.

The Project Community consists of all Contributors and Users of the
Project. Contributors work on behalf of and are responsible to the
larger Project Community and we strive to keep the barrier between
Contributors and Users as low as possible.

The Project is not a legal entity, nor does it currently have any formal
relationships with legal entities.


Governance
==========

This section describes the governance and leadership model of The
Project.

The foundations of Project governance are:

-  openness and transparency
-  active contribution
-  institutional neutrality


Traditionally, Project leadership was provided by a subset of Contributors,
called Core Developers, whose active and consistent contributions have been
recognized by their receiving “commit rights” to the Project GitHub
repositories. In general, all Project decisions are made through consensus among
the Core Developers with input from the Community.

While this approach has served us well, as the Project grows, we see a need for
a more formal governance model. The SciPy Core Developers expressed a
preference for a leadership model which includes a BDFL (Benevolent Dictator
for Life). Therefore, moving forward The Project leadership will consist of a
BDFL and Steering Council.

BDFL
----

The Project will have a BDFL (Benevolent Dictator for Life), who is currently
Pauli Virtanen. As Dictator, the BDFL has the authority to make all final
decisions for The Project. As Benevolent, the BDFL, in practice, chooses to
defer that authority to the consensus of the community discussion channels and
the Steering Council (see below). It is expected, and in the past has been the
case, that the BDFL will only rarely assert his/her final authority. Because
rarely used, we refer to BDFL’s final authority as a “special” or “overriding”
vote. When it does occur, the BDFL override typically happens in situations
where there is a deadlock in the Steering Council or if the Steering Council
asks the BDFL to make a decision on a specific matter. To ensure the
benevolence of the BDFL, The Project encourages others to fork the project if
they disagree with the overall direction the BDFL is taking. The BDFL may
delegate his/her authority on a particular decision or set of decisions to
any other Council member at his/her discretion.

The BDFL can appoint his/her successor, but it is expected that the Steering
Council would be consulted on this decision. If the BDFL is unable to appoint a
successor, the Steering Council will make this decision - preferably by
consensus, but if needed, by a majority vote.

Note that the BDFL can step down at any time, and acting in good faith, will
also listen to serious calls to do so. Also note that the BDFL is more a role
for fallback decision making rather than that of a director/CEO.

Steering Council
----------------

The Project will have a Steering Council that consists of Project Contributors
who have produced contributions that are substantial in quality and quantity,
and sustained over at least one year. The overall role of the Council is to
ensure, through working with the BDFL and taking input from the Community, the
long-term well-being of the project, both technically and as a community.

The Council will have a Chair, who is tasked with keeping the organizational
aspects of the functioning of the Council and the Project on track. The
Council will also appoint a Release Manager for the Project, who has final
responsibility for one or more releases.

During the everyday project activities, Council Members participate in all
discussions, code review, and other project activities as peers with all other
Contributors and the Community. In these everyday activities, Council Members
do not have any special power or privilege through their membership on the
Council. However, it is expected that because of the quality and quantity of
their contributions and their expert knowledge of the Project Software and
Services, Council Members will provide useful guidance, both technical and
in terms of project direction, to potentially less experienced contributors.

The Steering Council and its Members play a special role in certain situations.
In particular, the Council may:

-   Make decisions about the overall scope, vision, and direction of the
    project.
-   Make decisions about strategic collaborations with other organizations or
    individuals.
-   Make decisions about specific technical issues, features, bugs, and pull
    requests. They are the primary mechanism of guiding the code review process
    and merging pull requests.
-   Make decisions about the Services that are run by The Project and manage
    those Services for the benefit of the Project and Community.
-   Make decisions when regular community discussion does not produce consensus
    on an issue in a reasonable time frame.
-  Update policy documents, such as this one.

Council membership
~~~~~~~~~~~~~~~~~~

To become eligible for being a Steering Council Member, an individual must be a
Project Contributor who has produced contributions that are substantial in
quality and quantity, and sustained over at least one year. Potential Council
Members are nominated by existing Council members and voted upon by the
existing Council after asking if the potential Member is interested and willing
to serve in that capacity. The Council will be initially formed from the set of
existing Core Developers who, as of January 2017, have been significantly
active over the last two years.

When considering potential Members, the Council will look at candidates with a
comprehensive view of their contributions. This will include, but is not limited
to, code, code review, infrastructure work, mailing list and chat participation,
community help/building, education and outreach, design work, etc. We are
deliberately not setting arbitrary quantitative metrics (like “100 commits in
this repo”) to avoid encouraging behavior that plays to the metrics rather than
the project’s overall well-being. We want to encourage a diverse array of
backgrounds, viewpoints, and talents in our team, which is why we explicitly do
not define code as the sole metric on which council membership will be
evaluated.

If a Council Member becomes inactive in the project for a period of one year,
they will be considered for removal from the Council. Before removal, inactive
Member will be approached to see if they plan on returning to active
participation. If not, they will be removed immediately upon a Council
vote. If they plan on returning to active participation soon, they will be
given a grace period of one year. If they don’t return to active participation
within that time period they will be removed by vote of the Council without
further grace period. All former Council Members can be considered for
membership again at any time in the future, like any other Project Contributor.
Retired Council Members will be listed on the project website, acknowledging
the period during which they were active in the Council.

The Council reserves the right to eject current Members, other than the BDFL,
if they are deemed to be actively harmful to the project’s well-being, and
attempts at communication and conflict resolution have failed.

A list of current Steering Council Members is maintained at the
page :ref:`governance-people`.

Council Chair
~~~~~~~~~~~~~

The Chair will be appointed by the Steering Council. The Chair can stay on as
long as he/she wants, but may step down at any time and will listen to
serious calls to do so (similar to the BDFL role). The Chair will be
responsible for:

- Starting a review of the technical direction of the project (as captured by
  the :ref:`scipy-roadmap`) bi-yearly, around mid-April and mid-October.
- At the same times of the year, summarizing any relevant
  organizational updates and issues in the preceding period, and asking for
  feedback/suggestions on the mailing list.
- Ensuring the composition of the Steering Council stays current.
- Ensuring matters discussed in private by the Steering Council get
  summarized on the mailing list to keep the Community informed.
- Ensuring other important organizational documents (e.g., Code of Conduct,
  Fiscal Sponsorship Agreement) stay current after they are added.

Release Manager
~~~~~~~~~~~~~~~

The Release Manager has final responsibility for making a release.  This
includes:

- Proposing of and deciding on the timing of a release.
- Determining the content of a release in case there is no consensus on a
  particular change or feature.
- Creating the release and announcing it on the relevant public channels.

For more details on what those responsibilities look like in practice, see
:ref:`making-a-release`.

Conflict of interest
~~~~~~~~~~~~~~~~~~~~

It is expected that the BDFL and Council Members will be employed at a wide
range of companies, universities, and non-profit organizations. Because of this,
it is possible that Members will have a conflict of interest. Such conflicts of
interest include, but are not limited to:

-   Financial interest, such as investments, employment or contracting work,
    outside of The Project that may influence their work on The Project.
-   Access to proprietary information of their employer that could potentially
    leak into their work with the Project.

All members of the Council, BDFL included, shall disclose to the rest of the
Council any conflict of interest they may have. Members with a conflict of
interest in a particular issue may participate in Council discussions on that
issue, but must recuse themselves from voting on the issue. If the BDFL has
recused his/herself for a particular decision, the Council will appoint a
substitute BDFL for that decision.

Private communications of the Council
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unless specifically required, all Council discussions and activities will be
public and done in collaboration and discussion with the Project Contributors
and Community. The Council will have a private mailing list that will be used
sparingly and only when a specific matter requires privacy. When private
communications and decisions are needed, the Council will do its best to
summarize those to the Community after removing personal/private/sensitive
information that should not be posted to the public internet.

Council decision making
~~~~~~~~~~~~~~~~~~~~~~~

If it becomes necessary for the Steering Council to produce a formal
decision, then they will use a form of the `Apache Foundation voting
process <https://www.apache.org/foundation/voting.html>`_. This is a
formalized version of consensus, in which +1 votes indicate agreement,
-1 votes are vetoes (and must be accompanied with a rationale, as
above), and one can also vote fractionally (e.g. -0.5, +0.5) if one
wishes to express an opinion without registering a full veto. These
numeric votes are also often used informally as a way of getting a
general sense of people's feelings on some issue, and should not
normally be taken as formal votes. A formal vote only occurs if
explicitly declared, and if this does occur, then the vote should be held
open for long enough to give all interested Council Members a chance to
respond -- at least one week.

In practice, we anticipate that for most Steering Council decisions
(e.g., voting in new members) a more informal process will suffice.


Institutional Partners and funding
==================================

The Steering Council is the primary leadership for the project. No
outside institution, individual, or legal entity has the ability to own,
control, usurp, or influence the project other than by participating in
the Project as Contributors and Council Members. However, because
institutions can be an important funding mechanism for the project, it
is important to formally acknowledge institutional participation in the
project. These are Institutional Partners.

An Institutional Contributor is any individual Project Contributor who
contributes to the project as part of their official duties at an
Institutional Partner. Likewise, an Institutional Council Member is any
Project Steering Council Member who contributes to the project as part
of their official duties at an Institutional Partner.

With these definitions, an Institutional Partner is any recognized legal
entity in any country that employs at least 1 Institutional Contributor or
Institutional Council Member. Institutional Partners can be for-profit or
non-profit entities.

Institutions become eligible to become an Institutional Partner by
employing individuals who actively contribute to The Project as part of
their official duties. To state this another way, the only way for a
Partner to influence the project is by actively contributing to the open
development of the project, in equal terms to any other member of the
community of Contributors and Council Members. Merely using Project
Software in institutional context does not allow an entity to become an
Institutional Partner. Financial gifts do not enable an entity to become
an Institutional Partner. Once an institution becomes eligible for
Institutional Partnership, the Steering Council must nominate and
approve the Partnership.

If, at some point, an existing Institutional Partner stops having any
contributing employees, then a one year grace period commences. If, at
the end of this one-year period, they continue not to have any
contributing employees, then their Institutional Partnership will
lapse, and resuming it will require going through the normal process
for new Partnerships.

An Institutional Partner is free to pursue funding for their work on The
Project through any legal means. This could involve a non-profit
organization raising money from private foundations and donors or a
for-profit company building proprietary products and services that
leverage Project Software and Services. Funding acquired by
Institutional Partners to work on The Project is called Institutional
Funding. However, no funding obtained by an Institutional Partner can
override the Steering Council. If a Partner has funding to do SciPy work
and the Council decides to not pursue that work as a project, the
Partner is free to pursue it on their own. However, in this situation,
that part of the Partner’s work will not be under the SciPy umbrella and
cannot use the Project trademarks in any way that suggests a formal
relationship.

Institutional Partner benefits are:

-  acknowledgement on the SciPy website and in talks
-  ability to acknowledge their own funding sources on the SciPy
   website and in talks
-  ability to influence the project through the participation of their
   Council Member
-  invitation of the Council Members to SciPy Developer Meetings

A list of current Institutional Partners is maintained at the page
:ref:`governance-people`.


Document history
================

https://github.com/scipy/scipy/commits/master/doc/source/dev/governance/governance.rst

Acknowledgements
================

Substantial portions of this document were adapted from the
`Jupyter/IPython project's governance document
<https://github.com/jupyter/governance/blob/master/governance.md>`_ and
`NumPy's governance document
<https://github.com/numpy/numpy/blob/master/doc/source/dev/governance/governance.rst>`_.

License
=======

To the extent possible under law, the authors have waived all
copyright and related or neighboring rights to the SciPy project
governance document, as per the `CC-0 public domain dedication / license
<https://creativecommons.org/publicdomain/zero/1.0/>`_.

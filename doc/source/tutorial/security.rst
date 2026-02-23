.. _security:

Security
========

The SciPy maintenance team is security conscious and has
mechanisms for resolving identified security issues.
In general, however, core maintainers are not security experts
(their expertise lying elsewhere in scientific domains and programming).

We follow the :doc:`numpy:reference/security` guidelines
in addition to this document.

Data and Network Access
-----------------------

For working with untrusted data the :doc:`numpy:reference/security`
guidelines apply.

The majority of SciPy does not establish a network connection.
The sole exception to this is the :mod:`scipy.datasets` where
we rely on the `Pooch <https://www.fatiando.org/pooch/latest/>`_
dependency to fetch data and use sha256 codes to verify downloads.

Releases
--------

To mitigate the risk of supply chain attacks (as of the v1.18.0 release)
we rely on the `scipy-release <https://github.com/scipy/scipy-release>`_
repository which is modeled after the
`Numpy release <https://github.com/numpy/numpy-release>`_ process.

General security guidelines on the scipy-release repo:

- We require commit history and audit log are easy to inspect

- We require branch protection etc. and apply best practices

- We require all release artifacts to be built in the scipy-release repo

- We do not allow self-hosted runners on this repository

Disclosure
----------

To report a security vulnerability on the SciPy repo, please use
`Tidelift <https://tidelift.com/docs/security>`_.

SciPy is not designed to be exposed directly to untrusted users. A user
who can freely execute SciPy (or Python) functions must be considered
to have the same privileges as the process/Python interpreter.

If one can already execute Python code, there are far worse things one can do
than use all available CPU cycles, or provoke a symptom of a bug in code like
use-after-free or a segfault. Therefore, while such issues may be bugs, they
are not security issues.

Before reporting a security issue, please consider and describe the attack
vector in detail - and in particular whether that attack vector assumes being
able to freely execute SciPy functions.

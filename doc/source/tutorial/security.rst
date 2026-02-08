Security
========

In general the core maintainers of SciPy are scientific experts in their
domains and are not security experts. We are security conscious and have
mechanisms for resolving identified security issues.

We follow the guidelines under the NumPy security
guidelines :doc:`numpy:reference/security` in addition to this document.

Data and Network Access
-----------------------

For working with untrusted data the :doc:`numpy:reference/security`
guidelines apply.

The majority of SciPy does not establish a network connection.
The sole exception to this is the :mod:`scipy.datasets` where
we rely on the `Pooch <https://www.fatiando.org/pooch/latest/>`_
dependency to fetch data and use sha256 codes to verify the downloads.

Releases
--------

To mitigate the risk of supply chain attacks (as of the v1.18.0 release)
we rely on the `release <https://github.com/scipy/scipy-release>`_
repository which is modeled after the
`Numpy release <https://github.com/numpy/numpy-release>`_ process.

General security guidelines on the scipy-release repo:

- We require commit history and audit log are easy to inspect

- We require branch protection etc. and apply best practices

- We required all release artifacts to be built in the scipy-release repo

- We do not allow self-hosted runners on this repository

Disclosure
----------

To report vulnerabilities please review the guidelines here to determine
whether the vulnerability indeed needs to be reported on the SciPy repo and not
an upstream dependency. If the detected vulnerability requires remediation in
an upstream dependency we ask that you report the disclosure to the upstream
dependency.

To report a security vulnerability on this repo please use
`Tidelift <https://tidelift.com/docs/security>`_

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

Vendored Dependencies
---------------------

Exceptions to the upstream dependency policy would be if a vulnerability was in
a vendored dependency. Vendored dependencies are included
in `subprojects <https://github.com/scipy/scipy/tree/main/subprojects>`_

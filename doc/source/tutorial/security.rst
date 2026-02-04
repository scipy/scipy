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
guidelines apply. The majority of SciPy does not establish a network
connection.

The sole exception to this is the :mod"`scipy.datasets` where
we rely on the `Pooch <https://www.fatiando.org/pooch/latest/>`_
dependency to fetch data and use sha256 codes to verify the downloads.

Releases
--------

To mitigate the risk of supply chain attacks we rely on the
`release <https://github.com/scipy/scipy-release>`_ repository
which is modeled after the
`Numpy release <https://github.com/numpy/numpy-release>`_ process.

General security guidelines on the scipy-release repo:
* We require a linear history, so commit history is easy to inspect
* We require branch protection etc. and apply best practices
* We perform only wheels.yml runs, and optionally a security-related linter
action
* We required all release artifacts to be built inside this repository on
GitHub Actions runners
* We do not allow self-hosted runners on this repository
* We allow cross-compiling provided we are able to afford the billing and have
adequate maintainer time
* We perform verification of the test suite passing after cross compilation
either on the `scipy repo <https://github.com/scipy/scipy/tree/main>`_ or under
QEMU

Disclosure
----------

To report vulnerabilities please review the guidelines here to determine
whether the vulnerability indeed needs to be reported on the SciPy repo and not
an upstream dependency. If the detected vulnerability requires remediation in
an upstream dependency we ask that you report the disclosure to the upstream
dependency.

To report a security vulnerability on this repo please use
`TideLift <https://tidelift.com/docs/security>`_

Vendored Dependencies
---------------------

Exceptions to the upstream dependency policy would be if a vulnerability was in
a vendored dependency. Vendored dependencies are included
in `subprojects <https://github.com/scipy/scipy/tree/main/subprojects>`_

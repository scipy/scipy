.. _scipy-licensing:

Licensing
=========
SciPy is distributed under the `modified (3-clause) BSD license`_.  All code,
documentation and other files added to SciPy by contributors is licensed under
this license, unless another license is explicitly specified in the source
code.  Contributors keep the copyright for code they wrote and submit for
inclusion to SciPy.

Other licenses that are compatible with the modified BSD license that SciPy
uses are 2-clause BSD, MIT and PSF.  Incompatible licenses are GPL, Apache and
custom licenses that require attribution/citation or prohibit use for
commercial purposes.

PRs are often submitted with content copied or derived from unlicensed code or
code from a default license that is not compatible with SciPy's license.
For instance, code published on StackOverflow is covered by a CC-BY-SA license,
which is not compatible due to the share-alike clause.
These contributions cannot be accepted for inclusion in SciPy unless the
original code author is willing to (re)license their code under the
modified BSD (or compatible) license.  If the original author agrees,
add a comment saying so to the source files and forward the relevant
communication to the scipy-dev mailing list.

Another common occurrence is for code to be translated or derived from code in
R, Octave (both GPL-licensed) or a commercial application.  Such code also
cannot be included in SciPy.  Simply implementing functionality with the same
API as found in R/Octave/... is fine though, as long as the author doesn't look
at the original incompatibly-licensed source code.

.. _modified (3-clause) BSD license: https://opensource.org/licenses/BSD-3-Clause

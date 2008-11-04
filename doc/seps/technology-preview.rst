=============================
Technology Previews for SciPy
=============================

:Author: Jarrod Millman
:Contact: millman@berkeley.edu
:Date: 2008-11-03


Executive summary
=================

Technology preview code is new code incorporated into the trunk of
SciPy, but may also appear in a future SciPy release at our option.
Such code is considered production grade and well-tested, but which
also has no guarantees of a stable API to enable further improvements
based on community feedback.

In this document, I propose that we create a new top-level subpackage
in ``scipy`` called ``preview``.  The ``scipy.preview`` subpackage
will mirror ``scipy``.  This subpackage will be part of the official
releases.  Thus all code included in the technology preview will be
available to earlier testers and adopters, but clearly marked as under
very active development with no guareentee of a stable API.

scipy.preview
=============

The ``scipy.preview`` subpackage serves to main needs:

#. A place to include code that breaks current APIs of code in ``scipy``
#. A place to develop new packages, modules, classess, and functions
   that have got a stable API yet.

Consider the following scenarios.

For example, the new ``scipy.spatial`` subpackage is all new code that
is well-tested and fairly high quality.  Since we haven't had the time
to work out the ideal API and there are many possible use cases that we
may not have considered, there is a reluctance to freeze the API at this
time.  With the new ``scipy.preview`` subpackage, we can get more input
and reach a wider audience by releasing the current code in the next
release.

Another example would be if we decided that there is a function signature
that we want to change (e.g., percentileofscore http://codereview.appspot.com/7913).
In this case, we could put a deprecation warning in ``scipy.stats.percentileofscore``
and place the new function in ``scipy.preview.stats.percentileofscore``.

Requirements for inclusion
==========================

Getting code into ``scipy.preview`` would require a high barrier for entry.
We need to work out the details, but new packages, modules, classes, and
functions would need to be started in a branch or in another location.  The
authors would need to propose their code for inclusion into scipy.preview.
We would have to have some general agreement (or consensus) that we want
to include this functionality.  We would have to agree in general to the name
and would require it to follow our style guide and include a reasonable
amount of tests.

Requirements for promotion
==========================

For code in ``scipy.preview`` to be promoted requires ....

Advantages
==========

 * clear indication to user's what the future of the project holds
 * clear path for getting code into ``scipy``
 * lowers the barriers for API changes
 * clearly distinquishes what is most in flux
 * hopefully increase frequency of releases
 * gets new code out faster

running tests
=============

``scipy.test('full')`` doesn't run technology preview tests.  To run them would
require ``scipy.preview.test('full')``.

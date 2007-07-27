# req'd file for SciPy package (see DEVELOPERS.txt)
# Fred Mailhot
# 2006-06-13

raise ImportError(
"""ann has been moved to scikits. Please install
scikits.learn instead, and change your import to the following:

from scickits.learn.machine import ann.

For informations about scikits, see:
http://projects.scipy.org/scipy/scikits/""")

__all__ = ['mlp','srn','rbf']


major = 0
minor = 2
micro = 0
release_level = 'alpha'

from __cvs_version__ import cvs_version
cvs_minor = cvs_version[-3]
cvs_serial = cvs_version[-1]

scipy_version = '%(major)d.%(minor)d.%(micro)d-%(release_level)s'\
                '-%(cvs_minor)d.%(cvs_serial)d' % (locals ())

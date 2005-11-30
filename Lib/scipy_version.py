major = 0
minor = 4
micro = 3

scipy_version = '%(major)d.%(minor)d.%(micro)d' % (locals ())

import os
svn_version_file = os.path.join(os.path.dirname(__file__),
                                '__svn_version__.py')

if os.path.isfile(svn_version_file):
    import imp
    svn = imp.load_module('scipy.__svn_version__',
                          open(svn_version_file),
                          svn_version_file,
                          ('.py','U',1))
    scipy_version += '.'+svn.version


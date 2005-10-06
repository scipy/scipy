major = 0
minor = 2
micro = 0
release_level = 'beta'

from __svn_version__ import svn_version

special_version = '%(major)d.%(minor)d.%(micro)d_%(svn_version)s'\
                  % (locals ())

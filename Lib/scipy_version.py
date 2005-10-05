major = 0
minor = 4
micro = 2

try:
    from __svn_version__ import version as svn_revision
    scipy_version = '%(major)d.%(minor)d.%(micro)d_%(svn_revision)s'\
                    % (locals ())
except ImportError,msg:
    svn_revision = 0
    scipy_version = '%(major)d.%(minor)d.%(micro)d' % (locals ())


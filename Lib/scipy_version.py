#scipy_version = "0.2.0"

from scipy_distutils.misc_util import update_version
scipy_version = update_version(\
    version_template = \
    '%(major)d.%(minor)d.%(micro)d-%(release_level)s-%(serial)d',
    release_level = 'alpha', # alpha, beta, canditate, or final
                             # if final -> alpha then major = major + 1
    )

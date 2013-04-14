#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import os
from glob import glob

tmpl_empty_clapack_pyf = '''
python module clapack
  usercode void empty_module(void) {}
  interface
    subroutine empty_module()
      intent(c) empty_module
    end subroutine empty_module
  end interface
end python module clapack
'''


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info

    config = Configuration('lapack',parent_package,top_path)

    lapack_opt = get_info('lapack_opt',notfound_action=2)

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[])
                      if k == 'ATLAS_INFO']+[None])[0]
    if atlas_version:
        print(('ATLAS version: %s' % atlas_version))

    target_dir = ''

    # flapack:
    config.add_extension('flapack',
                         sources=['flapack.pyf.src'],
                         depends=[__file__,'flapack_*.pyf.src'],
                         extra_info=lapack_opt
                         )

    # clapack:
    def get_clapack_source(ext, build_dir):
        name = ext.name.split('.')[-1]
        assert name == 'clapack', repr(name)
        if atlas_version is None:
            target = os.path.join(build_dir,target_dir,'clapack.pyf')
            from distutils.dep_util import newer
            if newer(__file__,target):
                f = open(target,'w')
                f.write(tmpl_empty_clapack_pyf)
                f.close()
        else:
            target = ext.depends[0]
            assert os.path.basename(target) == 'clapack.pyf.src'
        return target

    config.add_extension('clapack',
                         sources=[get_clapack_source],
                         depends=['clapack.pyf.src'],
                         extra_info=lapack_opt
                         )

    # calc_lwork:
    config.add_extension('calc_lwork',
                         sources=['calc_lwork.f'],
                         extra_info=lapack_opt
                         )

    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())

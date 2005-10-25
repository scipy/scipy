# make source distribution and copy it (including CHANGES file) to wml directory

import os

OUTPUT_PATH = os.path.join(os.environ['HOME'], 'wml' , 'sites', 'www.geus.ch')

os.spawnlp(os.P_WAIT,
           'python',
           'python',
           'setup.py',
           'sdist')

execfile(os.path.join('Lib', 'pysparse_version.py'))
os.spawnlp(os.P_WAIT,
           'cp',
          'cp',
          os.path.join('dist', 'pysparse-%s.tar.gz' % version),
          os.path.join(OUTPUT_PATH, 'files'))

os.spawnlp(os.P_WAIT,
           'cp',
           'cp',
           'CHANGES',
           os.path.join(OUTPUT_PATH, 'files', 'CHANGES_pysparse.txt'))

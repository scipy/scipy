
#!/usr/bin/env python

from distutils.core import setup, Extension

MOD = 'NDI_Segmenter'
setup(name=MOD, ext_modules=[Extension(MOD, sources=['Segmenter_EXT.c', 'Segmenter_IMPL.c'])])


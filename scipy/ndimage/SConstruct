# Last Change: Wed Mar 05 09:00 PM 2008 J
from os.path import join

from numpy.distutils.misc_util import get_numpy_include_dirs
from numscons import GetNumpyEnvironment

env = GetNumpyEnvironment(ARGUMENTS)

env.AppendUnique(CPPPATH = get_numpy_include_dirs())
env.AppendUnique(CPPPATH = 'src')

ndimage_src = ["nd_image.c", "ni_filters.c", "ni_fourier.c", "ni_interpolation.c",
               "ni_measure.c", "ni_morphology.c", "ni_support.c"]
env.NumpyPythonExtension('_nd_image', source = [join('src', i) for i in ndimage_src])

segment_src = ['Segmenter_EXT.c', 'Segmenter_IMPL.c']
env.NumpyPythonExtension('_segment', source = [join('src', 'segment', i) 
                                               for i in segment_src])

register_src = ['Register_EXT.c', 'Register_IMPL.c']
env.NumpyPythonExtension('_register', source = [join('src', 'register', i) 
                                                for i in register_src])

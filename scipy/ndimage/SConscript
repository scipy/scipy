# Last Change: Fri Oct 10 03:00 PM 2008 J
from os.path import join

from numscons import GetNumpyEnvironment

env = GetNumpyEnvironment(ARGUMENTS)

env.AppendUnique(CPPPATH = 'src')

ndimage_src = ["nd_image.c", "ni_filters.c", "ni_fourier.c", "ni_interpolation.c",
               "ni_measure.c", "ni_morphology.c", "ni_support.c"]
env.NumpyPythonExtension('_nd_image', source = [join('src', i) for i in ndimage_src])

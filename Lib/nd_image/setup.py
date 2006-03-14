from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):

    config = Configuration('nd_image', parent_package, top_path,
                           package_path='Lib')

    config.add_extension("_nd_image",
        sources=["Src/nd_image.c","Src/ni_filters.c",
                 "Src/ni_fourier.c","Src/ni_interpolation.c",
                 "Src/ni_measure.c","Src/numcompat.c",
                 "Src/ni_morphology.c","Src/ni_support.c"],
        include_dirs=['Src'],
    )
    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())

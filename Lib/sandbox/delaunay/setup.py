from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):

    config = Configuration('delaunay', parent_package, top_path)

    config.add_extension("_delaunay",
        sources=["_delaunay.cpp", "VoronoiDiagramGenerator.cpp",
            "delaunay_utils.cpp", "natneighbors.cpp"],
        include_dirs=['.'],
    )

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())

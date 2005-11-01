from scipy.distutils.core import setup
from scipy.distutils.misc_util import get_path, Configuration

def configuration(parent_package='', parent_path=None):
    local_path = get_path(__name__)
    config = Configuration('delaunay', parent_package, parent_path)

    config.add_extension("_delaunay", 
        sources=["_delaunay.cpp", "VoronoiDiagramGenerator.cpp", 
            "delaunay_utils.cpp", "natneighbors.cpp"],
        include_dirs=[local_path],
    )

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())

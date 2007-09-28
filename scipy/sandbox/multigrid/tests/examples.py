import gzip
from scipy.io import mmread


def read_matrix(filename):
    filename = "sample_data/" + filename
    if filename.endswith(".gz"):
        fid = gzip.open(filename)
    else:
        fid = open(filename)

    return mmread(fid).tocsr()


mesh2d_laplacians = ['torus.mtx.gz','rocker_arm_surface.mtx.gz',
                     '336_triangle_A.mtx.gz','336_triangle_B.mtx.gz']


all_examples = mesh2d_laplacians

if __name__ == '__main__':
    print "All Available Examples Are Listed Below\n"
    for filename in all_examples:
        print filename
        print repr(read_matrix(filename))
        print "\n"


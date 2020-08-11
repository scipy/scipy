# distutils: language = c++

from .cy_declaration cimport hello_world

def cy_hello_world(double[::1] x):
    hello_world(&x[0], x.shape[0])

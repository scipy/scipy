# distutils: language = c++

from .cy_declaration cimport hello_world

def cy_hello_world():
    hello_world()
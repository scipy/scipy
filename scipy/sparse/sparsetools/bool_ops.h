#ifndef BOOL_OPS_H
#define BOOL_OPS_H

/*
 * Functions to handle arithmetic operations on NumPy Bool values.
 */

#include <numpy/arrayobject.h>

class npy_bool_wrapper {
    public:
        char value;
        
        /* operators */
        operator char() const {
            if(value != 0) {
                return 1;
            } else {
                return 0;
            }
        }
        npy_bool_wrapper& operator=(const npy_bool_wrapper& x) {
            value = x;
            return (*this);
        }
        npy_bool_wrapper operator+(const npy_bool_wrapper& x) {
            return x || value ? 1 : 0;
        }
        /* inplace operators */
        npy_bool_wrapper operator+=(const npy_bool_wrapper& x) {
            value = x || value ? 1 : 0;
            return (*this);
        }
        npy_bool_wrapper operator*=(const npy_bool_wrapper& x) {
            value = value && x ? 1 : 0;
            return (*this);
        }
        /* constructors */
        npy_bool_wrapper() { 
            value = 0; 
        }
        template <class T>
        npy_bool_wrapper(T x) {
            value = x ? 1 : 0;
        }
};

#endif

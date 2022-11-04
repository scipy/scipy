#ifndef BOOL_OPS_H
#define BOOL_OPS_H
/*
 * Functions to handle arithmetic operations on NumPy Bool values.
 */
#include <numpy/arrayobject.h>

/*
 * A compiler time (ct) assert macro from 
 * http://www.pixelbeat.org/programming/gcc/static_assert.html
 * This is used to assure that npy_bool_wrapper is the right size.
 */
#define ct_assert(e) extern char (*ct_assert(void)) [sizeof(char[1 - 2*!(e)])]

class npy_bool_wrapper {
    private:
        char value;

    public:
        /* operators */
        operator char() const {
            return value;
        }
        npy_bool_wrapper& operator=(const npy_bool_wrapper& x) {
            value = x.value;
            return (*this);
        }
        npy_bool_wrapper operator+(const npy_bool_wrapper& x) {
            return value || x.value;
        }
        /* inplace operators */
        npy_bool_wrapper operator+=(const npy_bool_wrapper& x) {
            value = (value || x.value);
            return (*this);
        }
        npy_bool_wrapper operator*=(const npy_bool_wrapper& x) {
            value = (value && x.value);
            return (*this);
        }
        /* constructors */
        npy_bool_wrapper() { 
            value = 0; 
        }
        template <class T>
        npy_bool_wrapper(T x) {
            value = (x != 0);
        }
};

ct_assert(sizeof(char) == sizeof(npy_bool_wrapper));

#endif


#include "nd_image.h"

PyArrayObject *
NA_NewArray(void *buffer, NumarrayType type, int ndim, ...)
{
        int i;
        maybelong shape[MAXDIM];
        va_list ap;
        va_start(ap, ndim);
        for(i=0; i<ndim; i++)
                shape[i] = va_arg(ap, int);  /* literals will still be ints */
        va_end(ap);
        return NA_vNewArray(buffer, type, ndim, shape);
}

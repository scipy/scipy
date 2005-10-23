#ifndef FORTRAN_H
#define FORTRAN_H

#ifdef NOF77UNDERSCORE
#define F77(s) s
#else
#define F77(s) s##_
#endif

#endif


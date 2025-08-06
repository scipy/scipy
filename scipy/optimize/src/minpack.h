#ifndef __MINPACK_H
#define __MINPACK_H


void CHKDER(const int,const int,double*,const double*,const double*,const int,double*,const double*,const int,double*);
void HYBRD(int(*)(int*,double*,double*,int*),const int,double*,double*,const double,const int,const int,const int,const double,double*,const int,const double,const int,int*,int*,double*,const int,double*,const int,double*,double*,double*,double*,double*);
void HYBRJ(int(*)(int*,double*,double*,double*,int*,int*),const int,double*,double*,double*,const int,const double,const int,double*,const int,const double,const int,int*,int*,int*,double*,const int,double*,double*,double*,double*,double*);
void LMDIF(int(*)(int*,int*,double*,double*,int*),const int,const int,double*,double*,const double,const double,const double,const int,const double,double*,const int,const double,const int,int*,int*,double*,const int,int*,double*,double*,double*,double*,double*);
void LMDER(int(*)(int*,int*,double*,double*,double*,int*,int*),const int,const int,double*,double*,double*,const int,const double,const double,const double,const int,double*,const int,const double,const int,int*,int*,int*,int*,double*,double*,double*,double*,double*);
void LMSTR(int(*)(int*,int*,double*,double*,double*,int*),const int,const int,double*,double*,double*,const int,const double,const double,const double,const int,double*,const int,const double,const int,int*,int*,int*,int*,double*,double*,double*,double*,double*);

#endif

#define NO_IMPORT_ARRAY
#include "ufunc_extras.h"

extern void PyUFunc_f_ff_As_d_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2];
	char *ip1=args[0], *op1=args[1], *op2=args[2];
	intp n=dimensions[0];
	double to1, to2;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2) {
		((IntFunc_d_dd *)func)((double)*(float *)ip1, &to1, &to2);
		*(float *)op1 = (float)to1;
		*(float *)op2 = (float)to2;
	}
}

extern void PyUFunc_d_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2];
	char *ip1=args[0], *op1=args[1], *op2=args[2];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2) {
		((IntFunc_d_dd *)func)(*(double *)ip1, (double *)op1, (double *)op2);
	}
}

extern void PyUFunc_F_FF_As_D_DD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2];
	char *ip1=args[0], *op1=args[1], *op2=args[2];
	intp n=dimensions[0];
        Py_complex from1;
	Py_complex to1, to2;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2) {
	        from1.real = ((float *)ip1)[0]; from1.imag = ((float *)ip1)[1];
		((IntFunc_D_DD *)func)(from1, &to1, &to2);
		((float *)op1)[0] = (float)to1.real;
		((float *)op1)[1] = (float)to1.imag;
		((float *)op2)[0] = (float)to2.real;
		((float *)op2)[1] = (float)to2.imag;
	}
}

extern void PyUFunc_D_DD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2];
	char *ip1=args[0], *op1=args[1], *op2=args[2];
	intp n=dimensions[0];
        Py_complex from1, to1, to2;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2) {
	        from1.real = ((double *)ip1)[0]; 
                from1.imag = ((double *)ip1)[1];
		((IntFunc_D_DD *)func)(from1, &to1, &to2);
		((double *)op1)[0] = to1.real;
		((double *)op1)[1] = to1.imag;
		((double *)op2)[0] = to2.real;
		((double *)op2)[1] = to2.imag;
        }

}

extern void PyUFunc_f_FF_As_d_DD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2];
	char *ip1=args[0], *op1=args[1], *op2=args[2];
	intp n=dimensions[0];
	Py_complex to1, to2;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2) {
		((IntFunc_d_DD *)func)((double)*(float *)ip1, &to1, &to2);
		((float *)op1)[0] = (float)to1.real;
		((float *)op1)[1] = (float)to1.imag;
		((float *)op2)[0] = (float)to2.real;
		((float *)op2)[1] = (float)to2.imag;
	}
}

extern void PyUFunc_d_DD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2];
	char *ip1=args[0], *op1=args[1], *op2=args[2];
	intp n=dimensions[0];
        Py_complex to1, to2;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2) {
		((IntFunc_d_DD *)func)(*(double *)ip1, &to1, &to2);
		((double *)op1)[0] = to1.real;
		((double *)op1)[1] = to1.imag;
		((double *)op2)[0] = to2.real;
		((double *)op2)[1] = to2.imag;
        }

}



extern void PyUFunc_f_ffff_As_d_dddd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2],os3=steps[3],os4=steps[4];
	char *ip1=args[0], *op1=args[1], *op2=args[2], *op3=args[3], *op4=args[4];
	intp n=dimensions[0];
	double to1, to2, to3, to4;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
		((IntFunc_d_dddd *)func)((double)*(float *)ip1, &to1, &to2, &to3, &to4);
		*(float *)op1 = (float)to1;
		*(float *)op2 = (float)to2;
		*(float *)op3 = (float)to3;
		*(float *)op4 = (float)to4;
	}
}


extern void PyUFunc_d_dddd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2],os3=steps[3],os4=steps[4];
	char *ip1=args[0], *op1=args[1], *op2=args[2], *op3=args[3], *op4=args[4];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
		((IntFunc_d_dddd *)func)(*(double *)ip1, (double *)op1, (double *)op2, (double *)op3, (double *)op4);
	}
}

extern void PyUFunc_F_FFFF_As_D_DDDD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2],os3=steps[3],os4=steps[4];
	char *ip1=args[0], *op1=args[1], *op2=args[2], *op3=args[3], *op4=args[4];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1, to2, to3, to4;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
	        from1.real = ((float *)ip1)[0]; from1.imag = ((float *)ip1)[1];
		((IntFunc_D_DDDD *)func)(from1, &to1, &to2, &to3, &to4);
		((float *)op1)[0] = (float)to1.real;
		((float *)op1)[1] = (float)to1.imag;
		((float *)op2)[0] = (float)to2.real;
		((float *)op2)[1] = (float)to2.imag;
		((float *)op3)[0] = (float)to3.real;
		((float *)op3)[1] = (float)to3.imag;
		((float *)op4)[0] = (float)to4.real;
		((float *)op4)[1] = (float)to4.imag;
	}
}

extern void PyUFunc_f_ffff_As_D_DDDD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2],os3=steps[3],os4=steps[4];
	char *ip1=args[0], *op1=args[1], *op2=args[2], *op3=args[3], *op4=args[4];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1, to2, to3, to4;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
	        from1.real = ((float *)ip1)[0]; from1.imag = 0.0;
		((IntFunc_D_DDDD *)func)(from1, &to1, &to2, &to3, &to4);
		((float *)op1)[0] = (float)to1.real;
		((float *)op2)[0] = (float)to2.real;
		((float *)op3)[0] = (float)to3.real;
		((float *)op4)[0] = (float)to4.real;
	}
}


extern void PyUFunc_d_dddd_As_D_DDDD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2],os3=steps[3],os4=steps[4];
	char *ip1=args[0], *op1=args[1], *op2=args[2], *op3=args[3], *op4=args[4];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1, to2, to3, to4;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
	        from1.real = ((double *)ip1)[0]; from1.imag = 0.0;
		((IntFunc_D_DDDD *)func)(from1, &to1, &to2, &to3, &to4);
		((double *)op1)[0] = (double)to1.real;
		((double *)op2)[0] = (double)to2.real;
		((double *)op3)[0] = (double)to3.real;
		((double *)op4)[0] = (double)to4.real;
	}
}

extern void PyUFunc_D_DDDD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2],os3=steps[3],os4=steps[4];
	char *ip1=args[0], *op1=args[1], *op2=args[2], *op3=args[3], *op4=args[4];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1, to2, to3, to4;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
	        from1.real = ((double *)ip1)[0]; from1.imag = ((double *)ip1)[1];
		((IntFunc_D_DDDD *)func)(from1, &to1, &to2, &to3, &to4);
		((double *)op1)[0] = (double)to1.real;
		((double *)op1)[1] = (double)to1.imag;
		((double *)op2)[0] = (double)to2.real;
		((double *)op2)[1] = (double)to2.imag;
		((double *)op3)[0] = (double)to3.real;
		((double *)op3)[1] = (double)to3.imag;
		((double *)op4)[0] = (double)to4.real;
		((double *)op4)[1] = (double)to4.imag;
	}
}


extern void PyUFunc_f_FFFF_As_d_DDDD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2],os3=steps[3],os4=steps[4];
	char *ip1=args[0], *op1=args[1], *op2=args[2], *op3=args[3], *op4=args[4];
	intp n=dimensions[0];
	Py_complex to1, to2, to3, to4;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
		((IntFunc_d_DDDD *)func)((double)(*((float *)ip1)), &to1, &to2, &to3, &to4);
		((float *)op1)[0] = (float)to1.real;
		((float *)op1)[1] = (float)to1.imag;
		((float *)op2)[0] = (float)to2.real;
		((float *)op2)[1] = (float)to2.imag;
		((float *)op3)[0] = (float)to3.real;
		((float *)op3)[1] = (float)to3.imag;
		((float *)op4)[0] = (float)to4.real;
		((float *)op4)[1] = (float)to4.imag;
	}
}

extern void PyUFunc_d_DDDD(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],os1=steps[1],os2=steps[2],os3=steps[3],os4=steps[4];
	char *ip1=args[0], *op1=args[1], *op2=args[2], *op3=args[3], *op4=args[4];
	intp n=dimensions[0];
	Py_complex to1, to2, to3, to4;
	
	for(i=0; i<n; i++, ip1+=is1, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
		((IntFunc_d_DDDD *)func)(*((double *)ip1), &to1, &to2, &to3, &to4);
		((double *)op1)[0] = (double)to1.real;
		((double *)op1)[1] = (double)to1.imag;
		((double *)op2)[0] = (double)to2.real;
		((double *)op2)[1] = (double)to2.imag;
		((double *)op3)[0] = (double)to3.real;
		((double *)op3)[1] = (double)to3.imag;
		((double *)op4)[0] = (double)to4.real;
		((double *)op4)[1] = (double)to4.imag;
	}
}

extern void PyUFunc_ff_ff_As_dd_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3];
	char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3];
	intp n=dimensions[0];
	double to1, to2;
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2) {
		((IntFunc_dd_dd *)func)((double)*(float *)ip1, (double)*(float *)ip2, &to1, &to2);
		*(float *)op1 = (float)to1;
		*(float *)op2 = (float)to2;
	}
}

extern void PyUFunc_dd_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3];
	char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2) {
		((IntFunc_dd_dd *)func)(*(double *)ip1, *(double *)ip2, (double *)op1, (double *)op2);
	}
}


extern void PyUFunc_ff_ffff_As_dd_dddd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4],os4=steps[5];
	char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4], *op4=args[5];
	intp n=dimensions[0];
	double to1, to2, to3, to4;
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
		((IntFunc_dd_dddd *)func)((double)*(float *)ip1, (double)*(float *)ip2, &to1, &to2, &to3, &to4);
		*(float *)op1 = (float)to1;
		*(float *)op2 = (float)to2;
		*(float *)op3 = (float)to3;
		*(float *)op4 = (float)to4;
	}
}

extern void PyUFunc_dd_dddd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os1=steps[2],os2=steps[3],os3=steps[4],os4=steps[5];
	char *ip1=args[0], *ip2=args[1], *op1=args[2], *op2=args[3], *op3=args[4], *op4=args[5];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op1+=os1, op2+=os2, op3+=os3, op4+=os4) {
		((IntFunc_dd_dddd *)func)(*(double *)ip1, *(double *)ip2, (double *)op1, (double *)op2, (double *)op3, (double *)op4);
	}
}


extern void PyUFunc_fff_f_As_ddd_d(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],os=steps[3];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os) {
		*(float *)op = (float)((DoubleFunc_ddd_d *)func)((double)*(float *)ip1, (double)*(float *)ip2, (double)*(float *)ip3);
	}
}

extern void PyUFunc_ddd_d(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],os=steps[3];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os) {
		*(double *)op = ((DoubleFunc_ddd_d *)func)(*(double *)ip1, *(double *)ip2, *(double *)ip3);
	}
}

extern void PyUFunc_fff_ff_As_ddd_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],os1=steps[3],os2=steps[4];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op1=args[3], *op2=args[4];
	intp n=dimensions[0];
	double to1, to2;
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op1+=os1, op2+=os2) {
		((IntFunc_ddd_dd *)func)((double)*(float *)ip1, (double)*(float *)ip2, (double)*(float *)ip3, &to1, &to2);
		*(float *)op1 = (float)to1;
		*(float *)op2 = (float)to2;
	}
}

extern void PyUFunc_ddd_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],os1=steps[3],os2=steps[4];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op1=args[3], *op2=args[4];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op1+=os1, op2+=os2) {
          ((IntFunc_ddd_dd *)func)(*(double *)ip1, *(double *)ip2, *(double *)ip3, (double *)op1, (double *)op2);
	}
}


extern void PyUFunc_ff_f_As_id_d(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os=steps[2];
	char *ip1=args[0], *ip2=args[1], *op=args[2];
	intp n=dimensions[0];

	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os) {
		*(float *)op = (float)((DoubleFunc_id_d *)func)((int)*(float *)ip1, (double)*(float *)ip2);
	}
}

extern void PyUFunc_dd_d_As_id_d(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os=steps[2];
	char *ip1=args[0], *ip2=args[1], *op=args[2];
	intp n=dimensions[0];

	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os) {
		*(double *)op = ((DoubleFunc_id_d *)func)((int)*(double *)ip1, *(double *)ip2);
	}
}


extern void PyUFunc_ff_f_As_dD_D(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os=steps[2];
	char *ip1=args[0], *ip2=args[1], *op=args[2];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1;

	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os) {
	        from1.real = *((float*)ip2); from1.imag = 0.0;
		to1 = ((CmplxFunc_dD_D *)func)((double)*(float *)ip1, from1);
		((float *)op)[0] = (float)to1.real;
	}
}

extern void PyUFunc_dd_d_As_dD_D(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os=steps[2];
	char *ip1=args[0], *ip2=args[1], *op=args[2];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1;

	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os) {
	        from1.real = *((double*)ip2); from1.imag = 0.0;
		to1 = ((CmplxFunc_dD_D *)func)(*(double *)ip1, from1);
		((double *)op)[0] = to1.real;
	}
}


extern void PyUFunc_fF_F_As_dD_D(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os=steps[2];
	char *ip1=args[0], *ip2=args[1], *op=args[2];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1;

	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os) {
	        from1.real = ((float*)ip2)[0]; from1.imag = ((float *)ip2)[1];
		to1 = ((CmplxFunc_dD_D *)func)((double)*(float *)ip1, from1);
		((float *)op)[0] = (float)to1.real;
		((float *)op)[1] = (float)to1.imag;
	}
}

extern void PyUFunc_dD_D(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],os=steps[2];
	char *ip1=args[0], *ip2=args[1], *op=args[2];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1;

	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, op+=os) {
	        from1.real = ((double*)ip2)[0]; from1.imag = ((double *)ip2)[1];
		to1 = ((CmplxFunc_dD_D *)func)(*(double *)ip1, from1);
		((double *)op)[0] = (double)to1.real;
		((double *)op)[1] = (double)to1.imag;
	}
}


extern void PyUFunc_ffF_F_As_ddD_D(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],os=steps[3];
	char *ip1=args[0],*ip2=args[1],*ip3=args[2],*op=args[3];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1;

	for(i=0; i<n; i++, ip1+=is1,ip2+=is2,ip3+=is3,op+=os) {
	        from1.real = ((float*)ip3)[0]; from1.imag = ((float *)ip3)[1];
		to1 = ((CmplxFunc_ddD_D *)func)((double)*(float *)ip1, (double)*(float *)ip2, from1);
		((float *)op)[0] = (float)to1.real;
		((float *)op)[1] = (float)to1.imag;
	}
}

extern void PyUFunc_ddD_D(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],os=steps[3];
	char *ip1=args[0],*ip2=args[1],*ip3=args[2],*op=args[3];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1;

	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os) {
	        from1.real = ((double*)ip3)[0]; from1.imag = ((double *)ip3)[1];
		to1 = ((CmplxFunc_ddD_D *)func)(*(double *)ip1, *(double *)ip2, from1);
		((double *)op)[0] = (double)to1.real;
		((double *)op)[1] = (double)to1.imag;
	}
}


extern void PyUFunc_fffF_F_As_dddD_D(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os=steps[4];
	char *ip1=args[0],*ip2=args[1],*ip3=args[2],*ip4=args[3],*op=args[4];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1;

	for(i=0; i<n; i++, ip1+=is1,ip2+=is2,ip3+=is3,ip4+=is4,op+=os) {
	        from1.real = ((float*)ip4)[0]; from1.imag = ((float *)ip4)[1];
		to1 = ((CmplxFunc_dddD_D *)func)((double)*(float *)ip1, (double)*(float *)ip2, (double)*(float *)ip3, from1);
		((float *)op)[0] = (float)to1.real;
		((float *)op)[1] = (float)to1.imag;
	}
}

extern void PyUFunc_dddD_D(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os=steps[4];
	char *ip1=args[0],*ip2=args[1],*ip3=args[2],*ip4=args[3],*op=args[4];
	intp n=dimensions[0];
	Py_complex from1;
	Py_complex to1;

	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op+=os) {
	        from1.real = ((double*)ip4)[0]; from1.imag = ((double *)ip4)[1];
		to1 = ((CmplxFunc_dddD_D *)func)(*(double *)ip1, *(double *)ip2, *(double *)ip3, from1);
		((double *)op)[0] = (double)to1.real;
		((double *)op)[1] = (double)to1.imag;
	}
}


extern void PyUFunc_fff_f_As_iid_d(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],os=steps[3];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os) {
		*(float *)op = (float)((DoubleFunc_iid_d *)func)((int)*(float *)ip1, (int)*(float *)ip2, (double)*(float *)ip3);
	}
}

extern void PyUFunc_ddd_d_As_iid_d(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],os=steps[3];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *op=args[3];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, op+=os) {
		*(double *)op = ((DoubleFunc_iid_d *)func)((int)*(double *)ip1, (int)*(double *)ip2, *(double *)ip3);
	}
}

extern void PyUFunc_ffff_f_As_dddd_d(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os=steps[4];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op=args[4];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op+=os) {
		*(float *)op = (float)((DoubleFunc_dddd_d *)func)((double)*(float *)ip1, (double)*(float *)ip2, (double)*(float *)ip3, (double)*(float *)ip4);
	}
}

extern void PyUFunc_dddd_d(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os=steps[4];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3],*op=args[4];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op+=os) {
		*(double *)op = ((DoubleFunc_dddd_d *)func)(*(double *)ip1, *(double *)ip2, *(double *)ip3, *(double *)ip4);
	}
}

extern void PyUFunc_ffff_ff_As_dddd_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4],os2=steps[5];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4], *op2=args[5];
	intp n=dimensions[0];
	double to1;
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op1+=os1, op2+=os2) {
		*(float *)op1 = (float)((DoubleFunc_dddd_dd *)func)((double)*(float *)ip1, (double)*(float *)ip2, (double)*(float *)ip3, (double)*(float *)ip4, &to1);
		*(float *)op2 = (float)to1;
	}
}

extern void PyUFunc_dddd_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4],os2=steps[5];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4], *op2=args[5];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op1+=os1, op2+=os2) {
		*(double *)op1 = ((DoubleFunc_dddd_dd *)func)(*(double *)ip1, *(double *)ip2, *(double *)ip3, *(double *)ip4, (double *)op2);
	}
}

extern void PyUFunc_fffff_ff_As_ddddd_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],is5=steps[4], os1=steps[5],os2=steps[6];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *ip5=args[4], *op1=args[5], *op2=args[6];
	intp n=dimensions[0];
	double to1, to2;
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, ip5+=is5, op1+=os1, op2+=os2) {
		((IntFunc_ddddd_dd *)func)((double)*(float *)ip1, (double)*(float *)ip2, (double)*(float *)ip3, (double)*(float *)ip4, (double)*(float *)ip5, &to1, &to2);
		*(float *)op1 = (float)to1;
		*(float *)op2 = (float)to2;
	}
}

extern void PyUFunc_ddddd_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],is5=steps[4],os1=steps[5],os2=steps[6];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *ip5=args[4], *op1=args[5], *op2=args[6];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, ip5+=is5, op1+=os1, op2+=os2) {
		((IntFunc_ddddd_dd *)func)(*(double *)ip1, *(double *)ip2, *(double *)ip3, *(double *)ip4, *(double *)ip5, (double *)op1, (double *)op2);
	}
}


extern void PyUFunc_ffff_ff_As_dddi_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4],os2=steps[5];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4], *op2=args[5];
	intp n=dimensions[0];
	double to1;
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op1+=os1, op2+=os2) {
		*(float *)op1 = (float)((DoubleFunc_dddi_dd *)func)((double)*(float *)ip1, (double)*(float *)ip2, (double)*(float *)ip3, (int)*(float *)ip4, &to1);
		*(float *)op2 = (float)to1;
	}
}

extern void PyUFunc_dddd_dd_As_dddi_dd(char **args, intp *dimensions, intp *steps, void *func) {
	int i, is1=steps[0],is2=steps[1],is3=steps[2],is4=steps[3],os1=steps[4],os2=steps[5];
	char *ip1=args[0], *ip2=args[1], *ip3=args[2], *ip4=args[3], *op1=args[4], *op2=args[5];
	intp n=dimensions[0];
	
	for(i=0; i<n; i++, ip1+=is1, ip2+=is2, ip3+=is3, ip4+=is4, op1+=os1, op2+=os2) {
		*(double *)op1 = ((DoubleFunc_dddi_dd *)func)(*(double *)ip1, *(double *)ip2, *(double *)ip3, (int)*(double *)ip4, (double *)op2);
	}
}














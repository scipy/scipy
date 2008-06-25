
extern "C" {
	double f1_pypy_g_bar(long l_n_1);
	double f2_pypy_g_bar(double l_n_5);
}


static long foo_1(long x)
{
	return f1_pypy_g_bar(x);
}

typedef long Func_1(long);

static void
PyUFunc_1(char **args, npy_intp *dimensions, npy_intp *steps, void *func)
{
	/* printf("PyUFunc_1\n"); */
	
	npy_intp n = dimensions[0];
	npy_intp is0 = steps[0];
	npy_intp os = steps[1];
	char *ip0 = args[0];
	char *op = args[1];
	Func_1 *f = (Func_1 *) func;
	npy_intp i;
	
	for(i = 0; i < n; i++, ip0 += is0, op += os) {
		long *in1 = (long *)ip0;
		long *out = (long *)op;
		
		*out = f(*in1);
	}
}

static double foo_2(double x)
{
	return f2_pypy_g_bar(x);
}

typedef double Func_2(double);

static void
PyUFunc_2(char **args, npy_intp *dimensions, npy_intp *steps, void *func)
{
	/* printf("PyUFunc_2\n"); */
	
	npy_intp n = dimensions[0];
	npy_intp is0 = steps[0];
	npy_intp os = steps[1];
	char *ip0 = args[0];
	char *op = args[1];
	Func_2 *f = (Func_2 *) func;
	npy_intp i;
	
	for(i = 0; i < n; i++, ip0 += is0, op += os) {
		double *in1 = (double *)ip0;
		double *out = (double *)op;
		
		*out = f(*in1);
	}
}


static PyUFuncGenericFunction foo_functions[] = {
	PyUFunc_1,
	PyUFunc_2,
};

static void *foo_data[] = {
	(void *) foo_1,
	(void *) foo_2,
};

static char foo_signatures[] = {
	NPY_LONG, NPY_LONG,   /* 1 */
	NPY_DOUBLE, NPY_DOUBLE,   /* 2 */
};

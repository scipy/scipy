#include<math.h>

typedef struct _complex cplx;

cplx cexp (cplx z)
{
	cplx res;
	double expre = exp(z.x);
	res.x = expre * cos(z.y);
	res.y = expre * sin(z.y);
	return res;
}

int main(void)
{
	cplx x = {10., 1.};
	cplx z;
	z = cexp(x);
	printf("%lf + i * %lf", z.x, z.y);
	return 0;
}
#include<math.h>

typedef struct _complex cplx;
typedef struct _cplxf{float x, y;} cplxf;

float cabsf(cplxf z)
{
	return sqrtf(z.x * z.x + z.y * z.y);
}

cplxf csqrtf(cplxf z)
{	
	cplxf ret;
	float r = cabsf(z);
	ret.x = sqrtf((r + z.x)/2.);
	ret.y = sqrtf((r - z.x)/2.);
	if (z.y < 0)
		ret.y = -ret.y;
	return ret;
}

cplx csqrt(cplx z)
{
	cplx ret;
	double r = cabs(z);
	ret.x = sqrt((r + z.x)/2.);
	ret.y = sqrt((r - z.x)/2.);
	if (z.y < 0)
		ret.y = -ret.y;
	return ret;
}
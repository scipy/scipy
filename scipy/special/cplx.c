#include<math.h>

typedef struct _complex cplx;

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

cplx cexp(cplx z)
{
	cplx res;
	double expre = exp(z.x);
	res.x = expre * cos(z.y);
	res.y = expre * sin(z.y);
	return res;
}

double __powidf2 (double x, int m)
{
  unsigned int n = m < 0 ? -m : m;
  double y = n % 2 ? x : 1;
  while (n >>= 1)
    {
      x = x * x;
      if (n % 2)
        y = y * x;
    }
  return m < 0 ? 1/y : y;
}

float __powisf2(float a, int b)
{
    const int recip = b < 0;
    float r = 1;
    while (1)
    {
        if (b & 1)
            r *= a;
        b /= 2;
        if (b == 0)
            break;
        a *= a;
    }
    return recip ? 1/r : r;
}

cplx clog(cplx z)
{
	cplx ret = {log(cabs(z)), atan2(z.y, z.x)};
	return ret;
}

cplx ccos(cplx z)
{
	cplx ret = {cos(z.x) * cosh(z.	y), - sin(z.x) * sinh(z.y)};
	return ret;
}

cplx csin(cplx z)
{
	cplx ret = {sin(z.x) * cosh(z.y), cos(z.x) * sinh(z.y)};
	return ret;
}

cplx cprod(cplx z1, cplx z2)
{
	cplx ret = {z1.x * z2.x - z1.y * z2.y, z1.x * z2.y + z2.x * z1.y};
	return ret;
}

cplx cpow (cplx base, cplx exponent)
{
	return cexp(cprod(exponent, clog(base))); 
}

double trunc(double x)
{
	if (x < 0)
		return -(int)abs(x);
	return (int)abs(x);
}


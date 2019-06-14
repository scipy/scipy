#ifdef ELLINT_POLY_REAL
#undef ELLINT_POLY_REAL
#endif
#define ELLINT_POLY_COMPLEX
#include "test_common.h"


#define NARGS	(3)


static int testwrapper(const EllInt_Num_t args[NARGS], double rerr,
		       EllInt_Num_t * restrict res)
{
    return cellint_RF(args[0], args[1], args[2], rerr, res);
}


int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const EllInt_Num_t eres[] = {
	1.3110287771461,
	1.8540746773014,
	1.8540746773014,
	0.79612586584234 - 1.2138566698365 * I,
	0.58408284167715,
	1.0441445654064,
	0.93912050218619 - 0.53296252018635 * I
    };
    const EllInt_Num_t args[][NARGS] = {
	{1.0, 2.0, 0.0},
	{I, -I, 0.0},
	{0.5, 1.0, 0.0},
	{-1.0 + I, I, 0.0},
	{2.0, 3.0, 4.0},
	{I, -I, 2.0},
	{-1.0 + I, I, 1.0 - I}
    };

    TESTARRAY_MAIN();
}

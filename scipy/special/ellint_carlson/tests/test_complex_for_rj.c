#ifdef ELLINT_POLY_REAL
#undef ELLINT_POLY_REAL
#endif
#define ELLINT_POLY_COMPLEX
#include "test_common.h"


#define NARGS	(4)


static int testwrapper(const EllInt_Num_t args[NARGS], double rerr,
		       EllInt_Num_t * restrict res)
{
    return cellint_RJ(args[0], args[1], args[2], args[3], rerr, res);
}


int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const EllInt_Num_t eres[] = {
	0.77688623778582,
	0.14297579667157,
	0.13613945827771 - 0.38207561624427 * I,
	1.6490011662711,
	0.94148358841220,
	1.8260115229009 + 1.2290661908643 * I,
	-0.61127970812028 - 1.0684038390007 * I,
	0.24723819703052,	/* Cauchy principal value */
	-0.12711230042964	/* Caucny principal value */
    };
    const EllInt_Num_t args[][NARGS] = {
	{0.0, 1.0, 2.0, 3.0},
	{2.0, 3.0, 4.0, 5.0},
	{2.0, 3.0, 4.0, -1.0 + I},
	{I, -I, 0.0, 2.0},
	{-1.0 + I, -1.0 - I, 1.0, 2.0},
	{I, -I, 0.0, 1.0 - I},
	{-1.0 + I, -1.0 - I, 1.0, -3.0 + I},
	{2.0, 3.0, 4.0, -0.5},	/* Cauchy principal value */
	{2.0, 3.0, 4.0, -5.0}	/* Cauchy principal value */
    };

    TESTARRAY_MAIN();
}

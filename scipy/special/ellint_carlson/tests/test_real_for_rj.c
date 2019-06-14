#ifdef ELLINT_POLY_COMPLEX
#undef ELLINT_POLY_COMPLEX
#endif
#define ELLINT_POLY_REAL
#include "test_common.h"


#define NARGS	(4)


static int testwrapper(const EllInt_Num_t args[NARGS], double rerr,
		       EllInt_Num_t * restrict res)
{
    return fellint_RJ(args[0], args[1], args[2], args[3], rerr, res);
}


int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const EllInt_Num_t eres[] = {
	0.77688623778582,
	0.14297579667157,
	0.24723819703052,	/* Cauchy principal value */
	-0.12711230042964	/* Caucny principal value */
    };
    const EllInt_Num_t args[][NARGS] = {
	{0.0, 1.0, 2.0, 3.0},
	{2.0, 3.0, 4.0, 5.0},
	{2.0, 3.0, 4.0, -0.5},	/* Cauchy principal value */
	{2.0, 3.0, 4.0, -5.0}	/* Cauchy principal value */
    };

    TESTARRAY_MAIN();
}

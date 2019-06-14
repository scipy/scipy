#ifdef ELLINT_POLY_COMPLEX
#undef ELLINT_POLY_COMPLEX
#endif
#define ELLINT_POLY_REAL
#include "test_common.h"


#define NARGS	(2)


static int testwrapper(const EllInt_Num_t args[NARGS], double rerr,
		       EllInt_Num_t * restrict res)
{
    return fellint_RC(args[0], args[1], rerr, res);
}


int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const EllInt_Num_t eres[] = {
	M_PI,
	M_LN2,
	M_LN2 / 3.0
    };
    const EllInt_Num_t args[][NARGS] = {
	{0.0, 0.25},
	{2.25, 2.0},
	{0.25, -2.0}
    };

    TESTARRAY_MAIN();
}

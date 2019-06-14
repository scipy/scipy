#ifdef ELLINT_POLY_COMPLEX
#undef ELLINT_POLY_COMPLEX
#endif
#define ELLINT_POLY_REAL
#include "test_common.h"


#define NARGS	(3)


static int testwrapper(const EllInt_Num_t args[NARGS], double rerr,
		       EllInt_Num_t * restrict res)
{
    return fellint_RD(args[0], args[1], args[2], rerr, res);
}


int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const EllInt_Num_t eres[] = {
	1.7972103521034,
	0.16510527294261
    };
    const EllInt_Num_t args[][NARGS] = {
	{0.0, 2.0, 1.0},
	{2.0, 3.0, 4.0}
    };

    TESTARRAY_MAIN();
}

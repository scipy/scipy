#ifdef ELLINT_POLY_COMPLEX
#undef ELLINT_POLY_COMPLEX
#endif
#define ELLINT_POLY_REAL
#include "test_common.h"


#define NARGS	(3)


static int testwrapper(const EllInt_Num_t args[NARGS], double rerr,
		       EllInt_Num_t * restrict res)
{
    return fellint_RG(args[0], args[1], args[2], rerr, res);
}


int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const EllInt_Num_t eres[] = {
	M_PI,
	1.7255030280692,
	1.0284758090288,
	0.5,
	0.0
    };
    const EllInt_Num_t args[][NARGS] = {
	{0.0, 16.0, 16.0},
	{2.0, 3.0, 4.0},
	{0.0, 0.0796, 4.0},
	{0.0, 0.0, 1.0},
	{0.0, 0.0, 0.0}
    };

    TESTARRAY_MAIN();
}

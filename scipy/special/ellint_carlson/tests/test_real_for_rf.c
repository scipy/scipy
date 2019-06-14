#ifdef ELLINT_POLY_COMPLEX
#undef ELLINT_POLY_COMPLEX
#endif
#define ELLINT_POLY_REAL
#include "test_common.h"


#define NARGS	(3)


static int testwrapper(const EllInt_Num_t args[NARGS], double rerr,
		       EllInt_Num_t * restrict res)
{
    return fellint_RF(args[0], args[1], args[2], rerr, res);
}


int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const EllInt_Num_t eres[] = {
	1.3110287771461,
	1.8540746773014,
	0.58408284167715
    };
    const EllInt_Num_t args[][NARGS] = {
	{1.0, 2.0, 0.0},
	{0.5, 1.0, 0.0},
	{2.0, 3.0, 4.0}
    };

    TESTARRAY_MAIN();
}

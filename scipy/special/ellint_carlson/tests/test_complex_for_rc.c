#ifdef ELLINT_POLY_REAL
#undef ELLINT_POLY_REAL
#endif
#define ELLINT_POLY_COMPLEX
#include "test_common.h"


#define NARGS	(2)


static int testwrapper(const EllInt_Num_t args[NARGS], double rerr,
                       EllInt_Num_t * restrict res)
{
    return cellint_RC(args[0], args[1], rerr, res);
}

int main(int argc, char * argv[])
{
    (void)argc;
    (void)argv;
    const EllInt_Num_t eres[] = {
	M_PI,
	M_LN2,
	1.1107207345396 - 1.1107207345396 * I,
	1.2260849569072 - 0.34471136988768 * I,
	M_LN2 / 3.0,
	0.77778596920447 + 0.19832484993429 * I
    };
    const EllInt_Num_t args[][NARGS] = {
	{0.0, 0.25},
	{2.25, 2.0},
	{0.0, I},
	{-I, I},
	{0.25, -2.0},
	{I, -1.0}
    };

    TESTARRAY_MAIN();
}

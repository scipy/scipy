#ifndef ELLINT_COMMON_HH_INCLUDED
#define ELLINT_COMMON_HH_INCLUDED


#include <cstddef>
#include <complex>


/* Namespace skeletons, to be populated in individual headers. */
namespace ellint_carlson
{
    /* Number of correction terms in polynomial evaluation. */
    static constexpr std::size_t NPOLY = 4u;

    /* Type-checking support. */
    namespace typing {}

    /* Arithmetic building blocks. */
    namespace arithmetic
    {
	namespace aux {}
    }

    /* Argument validation */
    namespace argcheck {}

    /* Configuration constants */
    namespace config
    {
	static constexpr unsigned int max_iter = 1000u;
	static constexpr double asym_zero_ul = 5.0e-14;
	static constexpr double asym_close_ul = 1.0e-9;
    }

    namespace constants
    {
	static constexpr long double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899863L;
	static constexpr long double ln4 = 1.38629436111989061883446424291635313615100026872051050824136L;

	static constexpr double RC_C[8] =
	{
	    80080.0,
	    0.0,
	    24024.0,	/* 3 / 10 */
	    11440.0,	/* 1 / 7 */
	    30030.0,	/* 3 / 8 */
	    32760.0,	/* 9 / 22 */
	    61215.0,	/* 159 / 208 */
	    90090.0	/* 9 / 8 */
	};

	static constexpr double RF_C1[4] =
	{
	    0.0,
	    -24024.0,
	    10010.0,
	    -5775.0
	};

	static constexpr double RF_C2[3] = {17160.0, -16380.0, 15015.0};

	static constexpr double RF_c33 = 6930.0;		/*  3 / 104 */

	static constexpr double RF_DENOM = 240240.0;

	static constexpr double RDJ_C1[4] = {0.0, -875160.0, 417690.0,
	                                     -255255.0};

	static constexpr double RDJ_C2[3] = {0.0, 680680.0, 306306.0};

	static constexpr double RDJ_C3[3] = {0.0, -706860.0, 675675.0};

	static constexpr double RDJ_C4[2] = {-556920.0, 612612.0};

	static constexpr double RDJ_C5[2] = {471240.0, -540540.0};

	static constexpr double RDJ_DENOM = 4084080.0;
    }  /* namespace constants */


    /* Exit statuses for the library functions. */
    enum class ExitStatus
    {
	success = 0,
	singular,
	underflow,
	overflow,
	n_iter,
	prec_loss,
	no_result,
	bad_args,
	bad_rerr,
	other,
	unused
    };

    inline static constexpr bool
    is_horrible(ExitStatus n)
    {
	return ( ( n == ExitStatus::no_result ) ||
	         ( n == ExitStatus::bad_args ) ||
		 ( n == ExitStatus::bad_rerr ) ||
		 ( n == ExitStatus::other ) );
    }

    inline static constexpr bool
    is_troublesome(ExitStatus n)
    {
	return ( !(( n == ExitStatus::success ) || is_horrible(n)) );
    }

    /* The upper-most namespace "ellint_carlson" will include the actual
     * functions (RC, RD, RF, etc.). */

    namespace util
    {
	/* Comparison function based on absolute value. */
	template<typename T>
	static inline bool abscmp(const T& a, const T& b)
	{
	    return ( std::abs(a) < std::abs(b) );
	}
    }
}


#endif /* ELLINT_COMMON_HH_INCLUDED */

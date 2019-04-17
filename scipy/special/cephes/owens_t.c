/* Copyright Benjamin Sobotta 2012
 *
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)
 */

/*
 * Reference:
 * Mike Patefield, David Tandy
 * FAST AND ACCURATE CALCULATION OF OWEN'S T-FUNCTION
 * Journal of Statistical Software, 5 (5), 1-25
 */
#include "mconf.h"
#include "_c99compat.h"

static const int SELECT_METHOD[] = {
    0, 0, 1, 12, 12, 12, 12, 12, 12, 12, 12, 15, 15, 15, 8,
    0, 1, 1, 2, 2, 4, 4, 13, 13, 14, 14, 15, 15, 15, 8,
    1, 1, 2, 2, 2, 4, 4, 14, 14, 14, 14, 15, 15, 15, 9,
    1, 1, 2, 4, 4, 4, 4, 6, 6, 15, 15, 15, 15, 15, 9,
    1, 2 , 2, 4, 4, 5 , 5, 7, 7, 16 ,16, 16, 11, 11, 10,
    1, 2 , 4, 4 , 4, 5 , 5, 7, 7, 16, 16, 16, 11, 11, 11,
    1, 2 , 3, 3, 5, 5 , 7, 7, 16, 16, 16, 16, 16, 11, 11,
    1, 2 , 3 , 3 , 5, 5, 17, 17, 17, 17, 16, 16, 16, 11, 11
};

static const double HRANGE[] = {0.02, 0.06, 0.09, 0.125, 0.26, 0.4, 0.6, 1.6,
    1.7, 2.33, 2.4, 3.36, 3.4, 4.8};

static const double ARANGE[] = {0.025, 0.09, 0.15, 0.36, 0.5, 0.9, 0.99999};

static const double ORD[] = {2, 3, 4, 5, 7, 10, 12, 18, 10, 20, 30, 0, 4, 7,
    8, 20, 0, 0};

static const int METHODS[] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4,
    5, 6};

static const double C[] = {
    0.99999999999999999999999729978162447266851932041876728736094298092917625009873,
    -0.99999999999999999999467056379678391810626533251885323416799874878563998732905968,
    0.99999999999999999824849349313270659391127814689133077036298754586814091034842536,
    -0.9999999999999997703859616213643405880166422891953033591551179153879839440241685,
    0.99999999999998394883415238173334565554173013941245103172035286759201504179038147,
    -0.9999999999993063616095509371081203145247992197457263066869044528823599399470977,
    0.9999999999797336340409464429599229870590160411238245275855903767652432017766116267,
    -0.999999999574958412069046680119051639753412378037565521359444170241346845522403274,
    0.9999999933226234193375324943920160947158239076786103108097456617750134812033362048,
    -0.9999999188923242461073033481053037468263536806742737922476636768006622772762168467,
    0.9999992195143483674402853783549420883055129680082932629160081128947764415749728967,
    -0.999993935137206712830997921913316971472227199741857386575097250553105958772041501,
    0.99996135597690552745362392866517133091672395614263398912807169603795088421057688716,
    -0.99979556366513946026406788969630293820987757758641211293079784585126692672425362469,
    0.999092789629617100153486251423850590051366661947344315423226082520411961968929483,
    -0.996593837411918202119308620432614600338157335862888580671450938858935084316004769854,
    0.98910017138386127038463510314625339359073956513420458166238478926511821146316469589567,
    -0.970078558040693314521331982203762771512160168582494513347846407314584943870399016019,
    0.92911438683263187495758525500033707204091967947532160289872782771388170647150321633673,
    -0.8542058695956156057286980736842905011429254735181323743367879525470479126968822863,
    0.73796526033030091233118357742803709382964420335559408722681794195743240930748630755,
    -0.58523469882837394570128599003785154144164680587615878645171632791404210655891158,
    0.415997776145676306165661663581868460503874205343014196580122174949645271353372263,
    -0.2588210875241943574388730510317252236407805082485246378222935376279663808416534365,
    0.1375535825163892648504646951500265585055789019410617565727090346559210218472356689,
    -0.0607952766325955730493900985022020434830339794955745989150270485056436844239206648,
    0.0216337683299871528059836483840390514275488679530797294557060229266785853764115,
    -0.00593405693455186729876995814181203900550014220428843483927218267309209471516256,
    0.0011743414818332946510474576182739210553333860106811865963485870668929503649964142,
    -1.489155613350368934073453260689881330166342484405529981510694514036264969925132E-4,
    9.072354320794357587710929507988814669454281514268844884841547607134260303118208E-6
};

static const double PTS[] = {
    0.35082039676451715489E-02, 0.31279042338030753740E-01,
    0.85266826283219451090E-01, 0.16245071730812277011E+00,
    0.25851196049125434828E+00, 0.36807553840697533536E+00,
    0.48501092905604697475E+00, 0.60277514152618576821E+00,
    0.71477884217753226516E+00, 0.81475510988760098605E+00,
    0.89711029755948965867E+00, 0.95723808085944261843E+00,
    0.99178832974629703586E+00
};

static const double WTS[] = {
    0.18831438115323502887E-01, 0.18567086243977649478E-01,
    0.18042093461223385584E-01, 0.17263829606398753364E-01,
    0.16243219975989856730E-01, 0.14994592034116704829E-01,
    0.13535474469662088392E-01, 0.11886351605820165233E-01,
    0.10070377242777431897E-01, 0.81130545742299586629E-02,
    0.60419009528470238773E-02, 0.38862217010742057883E-02,
    0.16793031084546090448E-02
};


static int get_method(double h, double a) {
    int ihint, iaint, i;

    ihint = 14;
    iaint = 7;

    for (i = 0; i < 14; i++) {
        if (h <= HRANGE[i]) {
            ihint = i;
            break;
        }
    }

    for (i = 0; i < 7; i++) {
        if (a <= ARANGE[i]) {
            iaint = i;
            break;
        }
    }
    return SELECT_METHOD[iaint * 15 + ihint];
}


static double owens_t_norm1(double x) {
    return erf(x / sqrt(2)) / 2;
}


static double owens_t_norm2(double x) {
    return erfc(x / sqrt(2)) / 2;
}


static double owensT1(double h, double a, double m) {
    int j = 1;
    int jj = 1;

    double hs = -0.5 * h * h;
    double dhs = exp(hs);
    double as = a * a;
    double aj = a / (2 * NPY_PI);
    double dj = expm1(hs);
    double gj = hs * dhs;

    double val = atan(a) / (2 * NPY_PI);

    while (1) {
	val += dj*aj / jj;

	if (m <= j) {
	    break;
	}
	j++;
	jj += 2;
	aj *= as;
	dj = gj - dj;
	gj *= hs / j;
    }

    return val;
}


static double owensT2(double h, double a, double ah, double m) {
    int i = 1;
    int maxi = 2 * m + 1;
    double hs = h * h;
    double as = -a * a;
    double y = 1.0 / hs;
    double val = 0.0;
    double vi = a*exp(-0.5 * ah * ah) / sqrt(2 * NPY_PI);
    double z = (ndtr(ah) - 0.5) / h;

    while (1) {
	val += z;
	if (maxi <= i) {
	    break;
	}
	z = y * (vi - i * z);
	vi *= as;
	i += 2;
    }
    val *= exp(-0.5 * hs) / sqrt(2 * NPY_PI);

    return val;
}


static double owensT3(double h, double a, double ah) {
    double aa, hh, y, vi, zi, result;
    int i;

    aa = a * a;
    hh = h * h;
    y = 1 / hh;

    vi = a * exp(-ah * ah/ 2) / sqrt(2 * NPY_PI);
    zi = owens_t_norm1(ah) / h;
    result = 0;

    for(i = 0; i<= 30; i++) {
        result += zi * C[i];
        zi = y * ((2 * i + 1) * zi - vi);
        vi *= aa;
    }

    result *= exp(-hh / 2) / sqrt(2 * NPY_PI);

    return result;
}


static double owensT4(double h, double a, double m) {
    double maxi, hh, naa, ai, yi, result;
    int i;

    maxi = 2 * m + 1;
    hh = h * h;
    naa = -a * a;

    i = 1;
    ai = a * exp(-hh * (1 - naa) / 2) / (2 * NPY_PI);
    yi = 1;
    result = 0;

    while (1) {
        result += ai * yi;

        if (maxi <= i) {
            break;
        }

        i += 2;
        yi = (1 - hh * yi) / i;
        ai *= naa;
    }

    return result;
}


static double owensT5(double h, double a) {
    double result, r, aa, nhh;
    int i;

    result = 0;
    r = 0;
    aa = a * a;
    nhh = -0.5 * h * h;

    for (i = 1; i < 14; i++) {
        r = 1 + aa * PTS[i - 1];
        result += WTS[i - 1] * exp(nhh * r) / r;
    }

    result *= a;

    return result;
}


static double owensT6(double h, double a) {
    double normh, y, r, result;

    normh = owens_t_norm2(h);
    y = 1 - a;
    r = atan2(y, (1 + a));
    result = normh * (1 - normh) / 2;

    if (r != 0) {
        result -= r * exp(-y * h * h / (2 * r)) / (2 * NPY_PI);
    }

    return result;
}


static double owens_t_dispatch(double h, double a, double ah) {
    int index, meth_code;
    double m, result;

    if (h == 0) {
        return atan(a) / (2 * NPY_PI);
    }
    if (a == 0) {
	return 0;
    }
    if (a == 1) {
        return owens_t_norm2(-h) * owens_t_norm2(h) / 2;
    }

    index = get_method(h, a);
    m = ORD[index];
    meth_code = METHODS[index];

    switch(meth_code) {
    case 1:
	result = owensT1(h, a, m);
	break;
    case 2:
	result = owensT2(h, a, ah, m);
	break;
    case 3:
	result = owensT3(h, a, ah);
	break;
    case 4:
	result = owensT4(h, a, m);
	break;
    case 5:
	result = owensT5(h, a);
	break;
    case 6:
	result = owensT6(h, a);
	break;
    default:
	result = NPY_NAN;
    }

    return result;
}


double owens_t(double h, double a) {
    double result, fabs_a, fabs_ah, normh, normah;

    if (cephes_isnan(h) || cephes_isnan(a)) {
        return NPY_NAN;
    }

    /* exploit that T(-h,a) == T(h,a) */
    h = fabs(h);

    /*
     * Use equation (2) in the paper to remap the arguments such that
     * h >= 0 and 0 <= a <= 1 for the call of the actual computation
     * routine.
     */
    fabs_a = fabs(a);
    fabs_ah = fabs_a * h;

    if (fabs_a == NPY_INFINITY) {
	/* See page 13 in the paper */
	result = owens_t_norm2(h);
    }
    else if (h == NPY_INFINITY) {
	result = 0;
    }
    else if (fabs_a <= 1) {
        result = owens_t_dispatch(h, fabs_a, fabs_ah);
    }
    else {
        if (fabs_ah <= 0.67) {
            normh = owens_t_norm1(h);
            normah = owens_t_norm1(fabs_ah);
            result = 0.25 - normh * normah -
                owens_t_dispatch(fabs_ah, (1 / fabs_a), h);
        }
        else {
            normh = owens_t_norm2(h);
            normah = owens_t_norm2(fabs_ah);
            result = (normh + normah) / 2 - normh * normah -
                owens_t_dispatch(fabs_ah, (1 / fabs_a), h);
        }
    }

    if (a < 0) {
	/* exploit that T(h,-a) == -T(h,a) */
	return -result;
    }

    return result;
}

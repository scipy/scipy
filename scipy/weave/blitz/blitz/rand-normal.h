/***************************************************************************
 * blitz/rand-normal.h    Random Gaussian (Normal) generator
 *
 * $Id: rand-normal.h 1414 2005-11-01 22:04:59Z cookedm $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This code was relicensed under the modified BSD license for use in SciPy
 * by Todd Veldhuizen (see LICENSE.txt in the weave directory).
 *
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************
 *
 * This generator transforms a (0,1] uniform distribution into
 * a Normal distribution.  Let u,v be (0,1] random variables. Then
 *
 *    x = sqrt(-2 ln v) cos(pi (2u-1))
 *
 * is N(0,1) distributed.
 *
 * Reference: Athanasios Papoulis, "Probability, random variables,
 * and stochastic processes," McGraw-Hill : Toronto, 1991.
 *
 ***************************************************************************/

#ifndef BZ_RAND_NORMAL_H
#define BZ_RAND_NORMAL_H

#ifndef BZ_RANDOM_H
 #include <blitz/random.h>
#endif

#ifndef BZ_RAND_UNIFORM_H
 #include <blitz/rand-uniform.h>
#endif

#include <math.h>

BZ_NAMESPACE(blitz)

template<typename P_uniform BZ_TEMPLATE_DEFAULT(Uniform)>
class Normal {

public:
    typedef double T_numtype;

    Normal(double mean = 0.0, double variance = 1.0, double = 0.0)
        : mean_(mean), sigma_(::sqrt(variance))
    { 
    }

    void randomize() 
    { 
        uniform_.randomize();
    }
  
    double random()
    { 
        double u, v;

        do {
            u = uniform_.random();
            v = uniform_.random();    
        } while (v == 0);

        return mean_ + sigma_ * ::sqrt(-2*::log(v)) * ::cos(M_PI * (2*u - 1));
    } 

private:
    double mean_, sigma_;
    P_uniform uniform_;
};

BZ_NAMESPACE_END

#endif // BZ_RAND_NORMAL_H


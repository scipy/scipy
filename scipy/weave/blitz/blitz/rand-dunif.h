/***************************************************************************
 * blitz/rand-dunif.h    Discrete uniform generator
 *
 * $Id: rand-dunif.h 1414 2005-11-01 22:04:59Z cookedm $
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
 ***************************************************************************/

#ifndef BZ_RAND_DUNIF_H
#define BZ_RAND_DUNIF_H

#ifndef BZ_RANDOM_H
 #include <blitz/random.h>
#endif

#ifndef BZ_RAND_UNIFORM_H
 #include <blitz/rand-uniform.h>
#endif

#include <math.h>

BZ_NAMESPACE(blitz)

template<typename P_uniform BZ_TEMPLATE_DEFAULT(Uniform)>
class DiscreteUniform {

public:
    typedef int T_numtype;
    typedef P_uniform T_uniform;

    DiscreteUniform(int low, int high, double=0)
        : low_(low), range_(high-low+1)
    { 
    }

    void randomize() 
    { 
        uniform_.randomize();
    }
  
    int random()
    { 
        return int(uniform_.random() * range_ + low_);
    } 

private:
    int low_, range_;
    T_uniform uniform_;
};

BZ_NAMESPACE_END

#endif // BZ_RAND_DUNIF_H


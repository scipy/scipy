
#ifndef CKDTREE_CPP_METHODS
#define CKDTREE_CPP_METHODS

#ifdef CKDTREE_METHODS_IMPL
#define CKDTREE_EXTERN extern "C"
#else
#define CKDTREE_EXTERN extern "C" 
struct ckdtree;
#endif 

extern npy_float64 infinity;
extern int number_of_processors;

#ifndef NPY_LIKELY
#define NPY_LIKELY(x) (x)
#endif

#ifndef NPY_UNLIKELY
#define NPY_UNLIKELY(x) (x)
#endif

#include <cmath>
#include <vector>
#include "ordered_pair.h"
#include "coo_entries.h"

#if defined(__GNUC__)

inline void 
prefetch_datapoint(const npy_float64 *x, const npy_intp m) 
{
    const int cache_line = 64;  // x86, amd64
    char *cur = (char*)x;
    char *end = (char*)(x+m);
    while (cur < end) { 
        __builtin_prefetch((void*)cur);
        cur += cache_line;
    }
}

#else

#if defined(_WIN32)

#include <xmmintrin.h>

inline void
prefetch_datapoint(const npy_float64 *x, const npy_intp m)
{
    const int cache_line = 64;  // x86, amd64
    char *cur = (char*)x;
    char *end = (char*)(x+m);
    while (cur < end) {
        _mm_prefetch((const char*)cur,_MM_HINT_T0);
        cur += cache_line;
    }
}

#else

#define prefetch_datapoint(x,y)

#endif // _WIN32
#endif // __GNUC__


/*
 * Utility functions
 * =================
 */
 
inline npy_float64 
dmax(const npy_float64 x, const npy_float64 y) 
{
    if (x > y)
        return x;
    else
        return y;
};

inline npy_float64 
dabs(const npy_float64 x)
{
    if (x > 0)
        return x;
    else
        return -x;
}

inline npy_float64 
wrap_distance(const npy_float64 x, const npy_float64 hb, const npy_float64 fb)
{
    npy_float64 x1;
    if (NPY_UNLIKELY(x < -hb)) x1 = fb + x;
    else if (NPY_UNLIKELY(x > hb)) x1 = x - fb;
    else x1 = x;
#if 0
    printf("dabs_b x : %g x1 %g\n", x, x1);
#endif
    return x1;
}

/*
 * Measuring distances
 * ===================
 */

#if !defined(__GNUC__)

inline npy_float64
sqeuclidean_distance_double(const npy_float64 *u, const npy_float64 *v, npy_intp n)
{
    npy_float64 s;
    npy_intp i;
    // manually unrolled loop, might be vectorized
    npy_float64 acc[4] = {0., 0., 0., 0.};
    for (i = 0; i < n/4; i += 4) {
        npy_float64 _u[4] = {u[i], u[i + 1], u[i + 2], u[i + 3]};
        npy_float64 _v[4] = {v[i], v[i + 1], v[i + 2], v[i + 3]};
        npy_float64 diff[4] = {_u[0] - _v[0], 
                               _u[1] - _v[1], 
                               _u[2] - _v[2],
                               _u[3] - _v[3]};
        acc[0] += diff[0] * diff[0];
        acc[1] += diff[1] * diff[1];
        acc[2] += diff[2] * diff[2];
        acc[3] += diff[3] * diff[3];
    }
    s = acc[0] + acc[1] + acc[2] + acc[3];
    if (i < n) {
        for(; i<n; ++i) {
            npy_float64 d = u[i] - v[i];
            s += d * d;
        }
    }
    return s;
} 

 
inline npy_float64 
_distance_p(const npy_float64 *x, const npy_float64 *y,
            const npy_float64 p, const npy_intp k,
            const npy_float64 upperbound)
{    
   /*
    * Compute the distance between x and y
    *
    * Computes the Minkowski p-distance to the power p between two points.
    * If the distance**p is larger than upperbound, then any number larger
    * than upperbound may be returned (the calculation is truncated).
    */
    
    npy_intp i;
    npy_float64 r;
    r = 0;
    if (NPY_LIKELY(p==2.)) {
        /*
        for (i=0; i<k; ++i) {
            z = x[i] - y[i];
            r += z*z;
            if (r>upperbound)
                return r;
        }*/
        return sqeuclidean_distance_double(x,y,k);
    } 
    else if (p==infinity) {
        for (i=0; i<k; ++i) {
            r = dmax(r,dabs(x[i]-y[i]));
            if (r>upperbound)
                return r;
        }
    } 
    else if (p==1.) {
        for (i=0; i<k; ++i) {
            r += dabs(x[i]-y[i]);
            if (r>upperbound)
                return r;
        }
    } 
    else {
        for (i=0; i<k; ++i) {
            r += std::pow(dabs(x[i]-y[i]),p);
            if (r>upperbound)
                 return r;
        }
    }
    return r;
} 


#else

/* Optimized path for GCC (and compatible compilers) using computed goto */
 
inline npy_float64 
_distance_p(const npy_float64 *x, const npy_float64 *y,
            const npy_float64 p, const npy_intp k,
            const npy_float64 upperbound)
{    
    
    void *jmptab[] = {&&minkowski, 
                      &&max_coordinate_difference, 
                      &&manhattan,
                      NULL,
                      &&euclid};
                      
    int dispatch = (p == infinity) | ((p == 1.)<<1) | ((p == 2.)<<2);
    goto *jmptab[dispatch];
    
euclid:
    {
        const npy_float64 *u = x;
        const npy_float64 *v = y;
        const npy_intp n = k;
        npy_float64 s;
        npy_intp i;
        npy_float64 acc[4] = {0., 0., 0., 0.};
        for (i = 0; i < n/4; i += 4) {
            npy_float64 _u[4] = {u[i], u[i + 1], u[i + 2], u[i + 3]};
            npy_float64 _v[4] = {v[i], v[i + 1], v[i + 2], v[i + 3]};
            npy_float64 diff[4] = {_u[0] - _v[0], 
                                   _u[1] - _v[1], 
                                   _u[2] - _v[2],
                                   _u[3] - _v[3]};
            acc[0] += diff[0] * diff[0];
            acc[1] += diff[1] * diff[1];
            acc[2] += diff[2] * diff[2];
            acc[3] += diff[3] * diff[3];
        }
        s = acc[0] + acc[1] + acc[2] + acc[3];
        if (i < n) {
            for(; i<n; ++i) {
                npy_float64 d = u[i] - v[i];
                s += d * d;
            }
        }
        return s;
    }
    
max_coordinate_difference:
    {
        npy_intp i;
        npy_float64 r = 0;
        for (i=0; i<k; ++i) {
            r = dmax(r,dabs(x[i]-y[i]));
            if (r>upperbound)
                return r;
        }
        return r;
    }

manhattan:
    {
        npy_intp i;
        npy_float64 r = 0;
        for (i=0; i<k; ++i) {
            r += dabs(x[i]-y[i]);
            if (r>upperbound)
                return r;
        }
        return r;
    }

minkowski:
    {
        npy_intp i;
        npy_float64 r = 0;
        for (i=0; i<k; ++i) {
            r += std::pow(dabs(x[i]-y[i]),p);
            if (r>upperbound)
                 return r;
        }
        return r;
    }
} 


#endif /* __GNUC__ */


inline npy_float64 
_distance_pp(const npy_float64 *x, const npy_float64 *y,
            const npy_float64 p, const npy_intp k,
            const npy_float64 upperbound, const ckdtreebox * box)
{    
   /*
    * Compute the distance between x and y
    *
    * Computes the Minkowski p-distance to the power p between two points.
   * If the distance**p is larger than upperbound, then any number larger
    * than upperbound may be returned (the calculation is truncated).
    *
    */
    
    npy_intp i;
    npy_float64 r, r1;
    r = 0;
    for (i=0; i<k; ++i) {
        r1 = wrap_distance(x[i] - y[i], box->hbox[i], box->fbox[i]);
        if (NPY_LIKELY(p==2.)) {
            r += r1 * r1;
        } else if (p==infinity) {
            r = dmax(r,r1);
        } else if (p==1.) {
            r += dabs(r1);
        } else {
            r += std::pow(dabs(r1),p);
        }
        if (r>upperbound) 
            return r;
    } 
    return r;
} 
 
static inline npy_float64 side_distance_from_min_max(
    const npy_float64 x,
    const npy_float64 min,
    const npy_float64 max,
    const npy_float64 p,
    const npy_float64 infinity,
    const npy_float64 hb,
    const npy_float64 fb
) {
    npy_float64 s, t, tmin, tmax;

    if (NPY_LIKELY(fb <= 0)) {
        /* non-periodic */
        s = 0; 
        t = x - max;
        if (t > s) {
            s = t;
        } else {
            t = min - x;
            if (t > s) s = t;
        }
    } else {
        /* periodic */
        s = 0;
        tmax = x - max;
        tmin = x - min;
        /* is the test point in this range */
        if(NPY_UNLIKELY(tmax >= 0 || tmin <= 0)) {
            /* no */
            tmax = dabs(tmax);
            tmin = dabs(tmin);
            /* tmin will be the closer edge */
            if(tmin > tmax) {
                t = tmin;
                tmin = tmax;
                tmax = t;
            }
            if(tmax < hb) {
                /* both edges are less than half a box. */
                /* no wrapping, use the closer edge */
                s = tmin;
            } else if(tmin > hb) {
                /* both edge are more than half a box. */
                /* wrapping on both edge, use the 
                 * wrapped further edge */
                s = fb - tmax;
            } else {
                /* the further side is wrapped */
                tmax = fb - tmax;
                if(tmin > tmax) {
                    s = tmax;
                } else {
                    s = tmin;
                }
            }
        } else {
            /* yes. min distance is 0 */
            s = 0; 
        }
    }
#if 0
    printf("sdfmm: s %g t %g x %g min %g max %g\n", s, t, x, min, max);
#endif 

    if (NPY_UNLIKELY(p == 1 || p == infinity)) {
        s = dabs(s);
    } else if (NPY_LIKELY(p == 2)) {
        s = s * s;
    } else {
        s = std::pow(s,p);
    }
        
    return s;
}

/* Build methods in C++ for better speed and GIL release */

CKDTREE_EXTERN PyObject*
build_ckdtree(ckdtree *self, npy_intp start_idx, npy_intp end_idx,
              npy_float64 *maxes, npy_float64 *mins, int _median, int _compact);


/* Query methods in C++ for better speed and GIL release */

CKDTREE_EXTERN PyObject*
query_knn(const ckdtree     *self, 
          npy_float64       *dd, 
          npy_intp          *ii, 
          const npy_float64 *xx,
          const npy_intp     n,
          const npy_intp     k, 
          const npy_float64  eps, 
          const npy_float64  p, 
          const npy_float64  distance_upper_bound,
          const ckdtreebox * box);
          
CKDTREE_EXTERN PyObject*
query_pairs(const ckdtree *self, 
            const npy_float64 r, 
            const npy_float64 p, 
            const npy_float64 eps,
            std::vector<ordered_pair> *results);
            
CKDTREE_EXTERN PyObject*
count_neighbors(const ckdtree *self,
                const ckdtree *other,
                npy_intp n_queries,
                npy_float64 *real_r,
                npy_intp *results,
                npy_intp *idx, 
                const npy_float64 p);
                               
CKDTREE_EXTERN PyObject*
query_ball_point(const ckdtree *self,
                 const npy_float64 *x,
                 const npy_float64 r,
                 const npy_float64 p,
                 const npy_float64 eps,
                 const npy_intp n_queries,
                 std::vector<npy_intp> **results);
                 
CKDTREE_EXTERN PyObject*                
query_ball_tree(const ckdtree *self,
                const ckdtree *other,
                const npy_float64 r,
                const npy_float64 p,
                const npy_float64 eps,
                std::vector<npy_intp> **results);                

CKDTREE_EXTERN PyObject*                 
sparse_distance_matrix(const ckdtree *self,
                       const ckdtree *other,
                       const npy_float64 p,
                       const npy_float64 max_distance,
                       std::vector<coo_entry> *results);

                  
#endif



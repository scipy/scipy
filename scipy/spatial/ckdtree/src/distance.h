#include "distance_base.h"

struct PlainDist1D {
    static inline const double side_distance_from_min_max(
        const ckdtree * tree, const double x,
        const double min,
        const double max,
        const ckdtree_intp_t k
        )
    {
        double s, t;
        s = 0;
        t = x - max;
        if (t > s) {
            s = t;
        } else {
            t = min - x;
            if (t > s) s = t;
        }
        return s;
    }
    static inline void
    interval_interval(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const ckdtree_intp_t k,
                        double *min, double *max)
    {
        /* Compute the minimum/maximum distance along dimension k between points in
         * two hyperrectangles.
         */
        *min = std::fmax(0., std::fmax(rect1.mins()[k] - rect2.maxes()[k],
                                       rect2.mins()[k] - rect1.maxes()[k]));

        *max =  std::fmax(rect1.maxes()[k] - rect2.mins()[k],
                          rect2.maxes()[k] - rect1.mins()[k]);
    }

    static inline double
    point_point(const ckdtree * tree,
               const double *x, const double *y,
                 const ckdtree_intp_t k) {
        return std::fabs(x[k] - y[k]);
    }
};

typedef BaseMinkowskiDistPp<PlainDist1D> MinkowskiDistPp;
typedef BaseMinkowskiDistPinf<PlainDist1D> MinkowskiDistPinf;
typedef BaseMinkowskiDistP1<PlainDist1D> MinkowskiDistP1;
typedef BaseMinkowskiDistP2<PlainDist1D> NonOptimizedMinkowskiDistP2;

/*
 * Measuring distances
 * ===================
 */


#ifndef CKDTREE_BLAS_DIST
#define CKDTREE_BLAS_DIST

double
sqeuclidean_distance_double_blas(const double * CKDTREE_RESTRICT u, 
                                 const double * CKDTREE_RESTRICT v, 
                                 const ckdtree_intp_t n);

#endif


#ifndef CKDTREE_NO_VECTORS // Vectorized version

inline static double
sqeuclidean_distance_double(const double * CKDTREE_RESTRICT u, 
                            const double * CKDTREE_RESTRICT v, 
                            const ckdtree_intp_t n)
{

    using vec2d  = double __attribute__ ((vector_size (2*sizeof(double))));

    // manually unrolled loop using GNU vector extensions
    const ckdtree_uintp_t un = static_cast<const ckdtree_uintp_t>(n);

    switch(un) {
        case 0: return 0.;
        case 1:
        {
            double d = u[0] - v[0];
            return d*d;
        }
        case 2:
        {
            vec2d _u = {u[0], u[1]};
            vec2d _v = {v[0], v[1]};
            vec2d diff = _u - _v;
            vec2d acc = diff * diff;
            return acc[0] + acc[1];
        }
        case 3:
        {
            vec2d _u[2] = {{u[0], u[1]}, {u[2], 0.0}};
            vec2d _v[2] = {{v[0], v[1]}, {v[2], 0.0}};
            vec2d diff[2] = {_u[0] - _v[0], _u[1] - _v[1]};
            vec2d acc = diff[0] * diff[0] + diff[1] * diff[1];
            return acc[0] + acc[1];
        }
        case 4:
        {
            vec2d _u[2] = {{u[0], u[1]}, {u[2], u[3]}};
            vec2d _v[2] = {{v[0], v[1]}, {v[2], v[3]}};
            vec2d diff[2] = {_u[0] - _v[0], _u[1] - _v[1]};
            vec2d acc = diff[0] * diff[0] + diff[1] * diff[1];
            return acc[0] + acc[1];
        }
        case 5:
        {
            vec2d _u[2] = {{u[0], u[1]}, {u[2], u[3]}};
            vec2d _v[2] = {{v[0], v[1]}, {v[2], v[3]}};
            vec2d diff[2] = {_u[0] - _v[0], _u[1] - _v[1]};
            vec2d acc = diff[0] * diff[0] + diff[1] * diff[1];
            auto diff4 = u[4] - v[4];
            return acc[0] + acc[1] + diff4 * diff4;
        }
        case 6:
        {
            vec2d d2[3];
            #pragma unroll
            for (int k = 0; k < 3; ++k) {
                vec2d _u = { u[2*k], u[2*k + 1] };
                vec2d _v = { v[2*k], v[2*k + 1] };
                auto diff = _u - _v;
                d2[k] = diff * diff;
            }
            vec2d acc = d2[0] + d2[1] + d2[2];
            return acc[0] + acc[1];
        }

        default: 
        {
            if (CKDTREE_LIKELY(n < 64)) {
                vec2d _u[2] = {{u[0], u[1]}, {u[2], u[3]}};
                vec2d _v[2] = {{v[0], v[1]}, {v[2], v[3]}};
                vec2d diff[2] = {_u[0] - _v[0], _u[1] - _v[1]};
                vec2d acc = diff[0] * diff[0] + diff[1] * diff[1];
                return acc[0] + acc[1] + sqeuclidean_distance_double(u+4, v+4, n-4);
            }
        }
    }
    return sqeuclidean_distance_double_blas(u, v, n);
}

#else // Scalar version when GNU vector extensions are not available



inline static double
sqeuclidean_distance_double(const double * CKDTREE_RESTRICT u, 
                            const double * CKDTREE_RESTRICT v, 
                            const ckdtree_intp_t n)
{

   const ckdtree_uintp_t un = static_cast<const ckdtree_uintp_t>(n);

    switch(un) {

        case 0: return 0.0;
        case 1:
        {
            const double d = u[0] - v[0];
            return d*d;
        }
        case 2:
        {
            const double d0 = u[0] - v[0];
            const double d1 = u[1] - v[1];
            return (d0*d0) + (d1*d1);
        }
        case 3:
        {
            const double diff[3] = {u[0] - v[0],
                                    u[1] - v[1],
                                    u[2] - v[2]};
            return (diff[0] * diff[0] +
                    diff[1] * diff[1] +
                    diff[2] * diff[2]);
        }
        case 4:
        {
            double d2[4];
            #pragma unroll
            for (int k = 0; k < 4; ++k) {
                auto diff = u[k] - v[k];
                d2[k] = diff * diff;
            }
            return (d2[0] + d2[1]) + (d2[2] + d2[3]);
        }
        case 5:
        {
            double d2[5];
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                auto diff = u[k] - v[k];
                d2[k] = diff * diff;
            }
            return (d2[0] + d2[1]) + (d2[2] + d2[3]) + d2[4];
        }
        case 6:
        {
            double d2[6];
            #pragma unroll
            for (int k = 0; k < 6; ++k) {
                auto diff = u[k] - v[k];
                d2[k] = diff * diff;
            }
            return (d2[0] + d2[1]) + (d2[2] + d2[3]) + (d2[4] + d2[5]);
        }
    
        default:
        {
            if (CKDTREE_LIKELY(n < 50)) {
                double d2[4];
                #pragma unroll
                for (int k = 0; k < 4; ++k) {
                    auto diff = u[k] - v[k];
                    d2[k] = diff * diff;
                }
                return (d2[0] + d2[1]) + (d2[2] + d2[3]) 
                    + sqeuclidean_distance_double(u+4, v+4, n-4);
            }
        }
    }
    return sqeuclidean_distance_double_blas(u, v, n);
}

#endif // SCALAR


struct MinkowskiDistP2: NonOptimizedMinkowskiDistP2 {
    static inline double
    point_point_p(const ckdtree * tree,
               const double *x, const double *y,
               const double p, const ckdtree_intp_t k,
               const double upperbound)
    {
        return sqeuclidean_distance_double(x, y, k);
    }
};

struct BoxDist1D {
    static inline void _interval_interval_1d (
        double min, double max,
        double *realmin, double *realmax,
        const double full, const double half
    )
    {
        /* Minimum and maximum distance of two intervals in a periodic box
         *
         * min and max is the nonperiodic distance between the near
         * and far edges.
         *
         * full and half are the box size and 0.5 * box size.
         *
         * value is returned in realmin and realmax;
         *
         * This function is copied from kdcount, and the convention
         * of is that
         *
         * min = rect1.min - rect2.max
         * max = rect1.max - rect2.min = - (rect2.min - rect1.max)
         *
         * We will fix the convention later.
         * */
        if (CKDTREE_UNLIKELY(full <= 0)) {
            /* A non-periodic dimension */
            /* \/     */
            if(max <= 0 || min >= 0) {
                /* do not pass though 0 */
                min = std::fabs(min);
                max = std::fabs(max);
                if(min < max) {
                    *realmin = min;
                    *realmax = max;
                } else {
                    *realmin = max;
                    *realmax = min;
                }
            } else {
                min = std::fabs(min);
                max = std::fabs(max);
                *realmax = std::fmax(max, min);
                *realmin = 0;
            }
            /* done with non-periodic dimension */
            return;
        }
        if(max <= 0 || min >= 0) {
            /* do not pass through 0 */
            min = std::fabs(min);
            max = std::fabs(max);
            if(min > max) {
                double t = min;
                min = max;
                max = t;
            }
            if(max < half) {
                /* all below half*/
                *realmin = min;
                *realmax = max;
            } else if(min > half) {
                /* all above half */
                *realmax = full - min;
                *realmin = full - max;
            } else {
                /* min below, max above */
                *realmax = half;
                *realmin = std::fmin(min, full - max);
            }
        } else {
            /* pass though 0 */
            min = -min;
            if(min > max) max = min;
            if(max > half) max = half;
            *realmax = max;
            *realmin = 0;
        }
    }
    static inline void
    interval_interval(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const ckdtree_intp_t k,
                        double *min, double *max)
    {
        /* Compute the minimum/maximum distance along dimension k between points in
         * two hyperrectangles.
         */
        _interval_interval_1d(rect1.mins()[k] - rect2.maxes()[k],
                    rect1.maxes()[k] - rect2.mins()[k], min, max,
                    tree->raw_boxsize_data[k], tree->raw_boxsize_data[k + rect1.m]);
    }

    static inline double
    point_point(const ckdtree * tree,
               const double *x, const double *y,
               const ckdtree_intp_t k)
    {
        double r1;
        r1 = wrap_distance(x[k] - y[k], tree->raw_boxsize_data[k + tree->m], tree->raw_boxsize_data[k]);
        r1 = std::fabs(r1);
        return r1;
    }

    static inline const double
    wrap_position(const double x, const double boxsize)
    {
        if (boxsize <= 0) return x;
        const double r = std::floor(x / boxsize);
        double x1 = x - r * boxsize;
        /* ensure result is within the box. */
        while(x1 >= boxsize) x1 -= boxsize;
        while(x1 < 0) x1 += boxsize;
        return x1;
    }

    static inline const double side_distance_from_min_max(
        const ckdtree * tree, const double x,
        const double min,
        const double max,
        const ckdtree_intp_t k
        )
    {
        double s, t, tmin, tmax;
        double fb = tree->raw_boxsize_data[k];
        double hb = tree->raw_boxsize_data[k + tree->m];

        if (fb <= 0) {
            /* non-periodic dimension */
            s = PlainDist1D::side_distance_from_min_max(tree, x, min, max, k);
            return s;
        }

        /* periodic */
        s = 0;
        tmax = x - max;
        tmin = x - min;
        /* is the test point in this range */
        if(CKDTREE_LIKELY(tmax < 0 && tmin > 0)) {
            /* yes. min distance is 0 */
            return 0;
        }

        /* no */
        tmax = std::fabs(tmax);
        tmin = std::fabs(tmin);

        /* make tmin the closer edge */
        if(tmin > tmax) { t = tmin; tmin = tmax; tmax = t; }

        /* both edges are less than half a box. */
        /* no wrapping, use the closer edge */
        if(tmax < hb) return tmin;

        /* both edge are more than half a box. */
        /* wrapping on both edge, use the
         * wrapped further edge */
        if(tmin > hb) return fb - tmax;

        /* the further side is wrapped */
        tmax = fb - tmax;
        if(tmin > tmax) return tmax;
        return tmin;
    }

    private:
    static inline double
    wrap_distance(const double x, const double hb, const double fb)
    {
        double x1;
        if (CKDTREE_UNLIKELY(x < -hb)) x1 = fb + x;
        else if (CKDTREE_UNLIKELY(x > hb)) x1 = x - fb;
        else x1 = x;
    #if 0
        printf("ckdtree_fabs_b x : %g x1 %g\n", x, x1);
    #endif
        return x1;
    }


};


typedef BaseMinkowskiDistPp<BoxDist1D> BoxMinkowskiDistPp;
typedef BaseMinkowskiDistPinf<BoxDist1D> BoxMinkowskiDistPinf;
typedef BaseMinkowskiDistP1<BoxDist1D> BoxMinkowskiDistP1;
typedef BaseMinkowskiDistP2<BoxDist1D> BoxMinkowskiDistP2;


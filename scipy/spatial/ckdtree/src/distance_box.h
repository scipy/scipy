struct BoxMinMaxDist1 {
    static inline void _interval_interval_1d (
        npy_float64 min, npy_float64 max,
        npy_float64 *realmin, npy_float64 *realmax,
        const npy_float64 full, const npy_float64 half
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
        if(max <= 0 || min >= 0) {
            /* do not pass through 0 */
            min = dabs(min);
            max = dabs(max);
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
                *realmin = dmin(min, full - max);
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
                        const npy_intp k,
                        npy_float64 *min, npy_float64 *max)
    {
        /* Compute the minimum/maximum distance along dimension k between points in
         * two hyperrectangles.
         */
        _interval_interval_1d(rect1.mins[k] - rect2.maxes[k],
                    rect1.maxes[k] - rect2.mins[k], min, max,
                    tree->raw_boxsize_data[k], tree->raw_boxsize_data[k + rect1.m]);
    }

    static inline npy_float64
    point_point(const ckdtree * tree, 
               const npy_float64 *x, const npy_float64 *y,
               const npy_intp k) 
    {
        npy_float64 r1;
        r1 = wrap_distance(x[k] - y[k], tree->raw_boxsize_data[k + tree->m], tree->raw_boxsize_data[k]);
        r1 = dabs(r1);
        return r1;
    }
};


typedef MinkowskiDistPp<BoxMinMaxDist1> BoxMinMaxDistPp;
typedef MinkowskiDistPinf<BoxMinMaxDist1> BoxMinMaxDistPinf;
typedef MinkowskiDistP1<BoxMinMaxDist1> BoxMinMaxDistP1;
typedef MinkowskiDistP2<BoxMinMaxDist1> BoxMinMaxDistP2;


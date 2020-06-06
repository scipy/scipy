template <typename Dist1D>
struct BaseMinkowskiDistPp {
    /* 1-d pieces
     * These should only be used if p != infinity
     */

    static inline void
    interval_interval_p(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const ckdtree_intp_t k, const double p,
                        double *min, double *max)
    {
        /* Compute the minimum/maximum distance along dimension k between points in
         * two hyperrectangles.
         */
        Dist1D::interval_interval(tree, rect1, rect2, k, min, max);
        *min = std::pow(*min, p);
        *max = std::pow(*max, p);
    }

    static inline void
    rect_rect_p(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const double p,
                        double *min, double *max)
    {
        *min = 0.;
        *max = 0.;
        for(ckdtree_intp_t i=0; i<rect1.m; ++i) {
            double min_, max_;

            Dist1D::interval_interval(tree, rect1, rect2, i, &min_, &max_);

            *min += std::pow(min_, p);
            *max += std::pow(max_, p);
        }
    }

    static inline double
    point_point_p(const ckdtree * tree,
               const double *x, const double *y,
               const double p, const ckdtree_intp_t k,
               const double upperbound)
    {
       /*
        * Compute the distance between x and y
        *
        * Computes the Minkowski p-distance to the power p between two points.
        * If the distance**p is larger than upperbound, then any number larger
        * than upperbound may be returned (the calculation is truncated).
        */

        ckdtree_intp_t i;
        double r, r1;
        r = 0;
        for (i=0; i<k; ++i) {
            r1 = Dist1D::point_point(tree, x, y, i);
            r += std::pow(r1, p);
            if (r>upperbound)
                 return r;
        }
        return r;
    }

    static inline double
    distance_p(const double s, const double p)
    {
        return std::pow(s,p);
    }
};

template <typename Dist1D>
struct BaseMinkowskiDistP1 : public BaseMinkowskiDistPp<Dist1D> {

    static inline void
    interval_interval_p(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const ckdtree_intp_t k, const double p,
                        double *min, double *max)
    {
        /* Compute the minimum/maximum distance along dimension k between points in
         * two hyperrectangles.
         */
        Dist1D::interval_interval(tree, rect1, rect2, k, min, max);
    }

    static inline void
    rect_rect_p(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const double p,
                        double *min, double *max)
    {
        *min = 0.;
        *max = 0.;
        for(ckdtree_intp_t i=0; i<rect1.m; ++i) {
            double min_, max_;

            Dist1D::interval_interval(tree, rect1, rect2, i, &min_, &max_);

            *min += min_;
            *max += max_;
        }
    }

    static inline double
    point_point_p(const ckdtree * tree,
               const double *x, const double *y,
               const double p, const ckdtree_intp_t k,
               const double upperbound)
    {
        ckdtree_intp_t i;
        double r;
        r = 0;
        for (i=0; i<k; ++i) {
            r += Dist1D::point_point(tree, x, y, i);
            if (r>upperbound)
                return r;
        }
        return r;
    }

    static inline double
    distance_p(const double s, const double p)
    {
        return s;
    }
};

template <typename Dist1D>
struct BaseMinkowskiDistPinf : public BaseMinkowskiDistPp<Dist1D> {
    static inline void
    interval_interval_p(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const ckdtree_intp_t k, double p,
                        double *min, double *max)
    {
        return rect_rect_p(tree, rect1, rect2, p, min, max);
    }

    static inline void
    rect_rect_p(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const double p,
                        double *min, double *max)
    {
        *min = 0.;
        *max = 0.;
        for(ckdtree_intp_t i=0; i<rect1.m; ++i) {
            double min_, max_;

            Dist1D::interval_interval(tree, rect1, rect2, i, &min_, &max_);

            *min = ckdtree_fmax(*min, min_);
            *max = ckdtree_fmax(*max, max_);
        }
    }

    static inline double
    point_point_p(const ckdtree * tree,
               const double *x, const double *y,
               const double p, const ckdtree_intp_t k,
               const double upperbound)
    {
        ckdtree_intp_t i;
        double r;
        r = 0;
        for (i=0; i<k; ++i) {
            r = ckdtree_fmax(r,Dist1D::point_point(tree, x, y, i));
            if (r>upperbound)
                return r;
        }
        return r;
    }
    static inline double
    distance_p(const double s, const double p)
    {
        return s;
    }
};

template <typename Dist1D>
struct BaseMinkowskiDistP2 : public BaseMinkowskiDistPp<Dist1D> {
    static inline void
    interval_interval_p(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const ckdtree_intp_t k, const double p,
                        double *min, double *max)
    {
        /* Compute the minimum/maximum distance along dimension k between points in
         * two hyperrectangles.
         */
        Dist1D::interval_interval(tree, rect1, rect2, k, min, max);
        *min *= *min;
        *max *= *max;
    }

    static inline void
    rect_rect_p(const ckdtree * tree,
                        const Rectangle& rect1, const Rectangle& rect2,
                        const double p,
                        double *min, double *max)
    {
        *min = 0.;
        *max = 0.;
        for(ckdtree_intp_t i=0; i<rect1.m; ++i) {
            double min_, max_;

            Dist1D::interval_interval(tree, rect1, rect2, i, &min_, &max_);
            min_ *= min_;
            max_ *= max_;

            *min += min_;
            *max += max_;
        }
    }
    static inline double

    point_point_p(const ckdtree * tree,
               const double *x, const double *y,
               const double p, const ckdtree_intp_t k,
               const double upperbound)
    {
        ckdtree_intp_t i;
        double r;
        r = 0;
        for (i=0; i<k; ++i) {
            double r1 = Dist1D::point_point(tree, x, y, i);
            r += r1 * r1;
            if (r>upperbound)
                return r;
        }
        return r;
    }
    static inline double
    distance_p(const double s, const double p)
    {
        return s * s;
    }
};


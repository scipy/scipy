#include <stdlib.h>
#include <stdio.h>

#include <fftw3.h>

enum type {
        DCT_I = 1,
        DCT_II = 2,
        DCT_III = 3,
        DCT_IV = 4,
};

int gen(int type, int sz)
{
        double *a, *b;
        fftw_plan p;
        int i, tp;

        a = fftw_malloc(sizeof(*a) * sz);
        if (a == NULL) {
                fprintf(stderr, "failure\n");
                exit(EXIT_FAILURE);
        }
        b = fftw_malloc(sizeof(*b) * sz);
        if (b == NULL) {
                fprintf(stderr, "failure\n");
                exit(EXIT_FAILURE);
        }

        for(i=0; i < sz; ++i) {
                a[i] = i;
        }

        switch(type) {
                case DCT_I:
                        tp = FFTW_REDFT00;
                        break;
                case DCT_II:
                        tp = FFTW_REDFT10;
                        break;
                case DCT_III:
                        tp = FFTW_REDFT01;
                        break;
                case DCT_IV:
                        tp = FFTW_REDFT11;
                        break;
                default:
                        fprintf(stderr, "unknown type\n");
                        exit(EXIT_FAILURE);
        }

        p = fftw_plan_r2r_1d(sz, a, b, tp, FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);

        for(i=0; i < sz; ++i) {
                printf("%f\n", b[i]);
        }
        fftw_free(b);
        fftw_free(a);

        return 0;
}

int main(int argc, char* argv[])
{
        int n, tp;

        if (argc < 3) {
                fprintf(stderr, "missing argument: program type n\n");
                exit(EXIT_FAILURE);
        }
        tp = atoi(argv[1]);
        n = atoi(argv[2]);

        gen(tp, n);
        fftw_cleanup();

        return 0;
}

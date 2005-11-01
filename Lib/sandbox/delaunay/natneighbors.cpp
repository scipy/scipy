
#include "delaunay_utils.h"
#include "natneighbors.h"

#include <stack>
#include <list>
#include <vector>
#include <math.h>

using namespace std;

NaturalNeighbors::NaturalNeighbors(int npoints, int ntriangles, double *x, double *y,
    double *centers, int *nodes, int *neighbors)
{
    this->npoints = npoints;
    this->ntriangles = ntriangles;
    this->x = x;
    this->y = y;
    this->centers = centers;
    this->nodes = nodes;
    this->neighbors = neighbors;

    this->radii2 = new double[ntriangles];
    for (int i=0; i<ntriangles; i++) {
        double x2 = x[INDEX3(nodes,i,0)] - INDEX2(centers,i,0);
        x2 = x2*x2;
        double y2 = y[INDEX3(nodes,i,0)] - INDEX2(centers,i,1);
        y2 = y2*y2;
        this->radii2[i] = x2 + y2;
    }

    this->rng_seed = 1234567890;
}

NaturalNeighbors::~NaturalNeighbors()
{
    delete[] this->radii2;
}

// Dumb implementation of the Park-Miller minimal standard PRNG
// I hate doing this, but I want complete control without > 5 lines of code.
double NaturalNeighbors::ranf()
{
    const double a=16807;
    const double m=2147483647.0;
    this->rng_seed = (long)(fmod(this->rng_seed * a, m));
    return this->rng_seed / m;
}

int NaturalNeighbors::find_containing_triangle(double targetx, double targety, int start_triangle)
{
    int final_triangle;
    final_triangle = walking_triangles(start_triangle, targetx, targety, 
        x, y, nodes, neighbors);
    return final_triangle;
}

double NaturalNeighbors::interpolate_one(double *z, double targetx, double targety,
    double defvalue, int &start_triangle)
{
    int t = find_containing_triangle(targetx, targety, start_triangle);
    if (t == -1) return defvalue;

    start_triangle = t;
    list<int> circumtri;

    circumtri.push_back(t);
    stack<int> stackA;
    stack<int> stackB;
    int tnew, i;

    for (i=0; i<3; i++) {
        tnew = INDEX3(this->neighbors, t, i);
        if (tnew != -1) {
            stackA.push(tnew);
            stackB.push(t);
        }
    }
    while (!stackA.empty()) {
        tnew = stackA.top();
        stackA.pop();
        t = stackB.top();
        stackB.pop();
        double d2 = (SQ(targetx - INDEX2(this->centers,tnew,0))
                   + SQ(targety - INDEX2(this->centers,tnew,1)));
        if (d2 <= this->radii2[tnew]) {
            // tnew is a circumtriangle of the target
            circumtri.push_back(tnew);
            for (i=0; i<3; i++) {
                int ti = INDEX3(this->neighbors, tnew, i);
                if ((ti != -1) && (ti != t)) {
                    stackA.push(ti);
                    stackB.push(tnew);
                }
            }
        }
    }

    list<int>::iterator it;
    double f = 0.0;
    double A = 0.0;

    for (it = circumtri.begin(); it != circumtri.end(); it++) {
        int t = *it;
        double vx = INDEX2(this->centers, t, 0);
        double vy = INDEX2(this->centers, t, 1);
        vector<double> c(6);
        for (int i=0; i<3; i++) {
            int j = EDGE0(i);
            int k = EDGE1(i);

            if (!circumcenter(
                    this->x[INDEX3(this->nodes, t, j)],
                    this->y[INDEX3(this->nodes, t, j)],
                    this->x[INDEX3(this->nodes, t, k)],
                    this->y[INDEX3(this->nodes, t, k)],
                    targetx, targety,
                    INDEX2(c, i, 0), INDEX2(c, i, 1))) {

                // bail out with the appropriate values if we're actually on a 
                // node
                if ((targetx == this->x[INDEX3(this->nodes, t, j)]) &&
                    (targety == this->y[INDEX3(this->nodes, t, j)])) {
                    return z[INDEX3(this->nodes, t, j)];
                } else if ((targetx == this->x[INDEX3(this->nodes, t, k)]) &&
                           (targety == this->y[INDEX3(this->nodes, t, k)])) {
                    return z[INDEX3(this->nodes, t, k)];
                } else {
                    // We're somewhere on one of the Delaunay lines
                    // Perturb the result a little bit and redo the calculation
                    double dx, dy, f;
                    dx = PERTURB_EPS * (this->y[INDEX3(this->nodes, t, j)] 
                                      - this->y[INDEX3(this->nodes, t, k)]);
                    dy = PERTURB_EPS * (this->x[INDEX3(this->nodes, t, j)] 
                                      - this->x[INDEX3(this->nodes, t, k)]);
                    f = interpolate_one(z, targetx+dx*ranf(), targety+dy*ranf(), 
                        defvalue, start_triangle);
                    return f;
                }
            }
        }
        for (int i=0; i<3; i++) {
            int j = EDGE0(i);
            int k = EDGE1(i);

            int q = INDEX3(this->nodes, t, i);
            double ati = signed_area(vx, vy, 
                                     INDEX2(c, j, 0), INDEX2(c, j, 1),
                                     INDEX2(c, k, 0), INDEX2(c, k, 1));
            A += ati;
            f += ati*z[q];
        }
    }
    f /= A;

    return f;
}

void NaturalNeighbors::interpolate_grid(double *z, 
    double x0, double x1, int xsteps,
    double y0, double y1, int ysteps,
    double *output,
    double defvalue, int start_triangle)
{
    int i, ix, iy, rowtri, coltri, tri;
    double dx, dy, targetx, targety;

    dx = (x1 - x0) / (xsteps-1);
    dy = (y1 - y0) / (ysteps-1);

    rowtri = 0;
    i = 0;
    for (iy=0; iy<ysteps; iy++) {
        targety = y0 + dy*iy;
        rowtri = find_containing_triangle(x0, targety, rowtri);
        tri = rowtri;
        for (ix=0; ix<xsteps; ix++) {
            targetx = x0 + dx*ix;
            coltri = tri;
            INDEXN(output, xsteps, iy, ix) = interpolate_one(z, targetx, targety,
                defvalue, coltri);
            if (coltri != -1) tri = coltri;
        }
    }
}

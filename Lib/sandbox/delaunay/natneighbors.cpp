
#include "delaunay_utils.h"
#include "natneighbors.h"

#include <stack>
#include <list>
#include <vector>
#include <set>
#include <math.h>
#include <iostream>
#include <iterator>

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
    bool test = false; //((targetx == 1.0));
    int t = find_containing_triangle(targetx, targety, start_triangle);
    if (t == -1) return defvalue;

    start_triangle = t;
    vector<int> circumtri;

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
        if (d2 < this->radii2[tnew]) {
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

    vector<int>::iterator it;
    double f = 0.0;
    double A = 0.0;
    double tA=0.0, yA=0.0, cA=0.0; // Kahan summation temps for A
    double tf=0.0, yf=0.0, cf=0.0; // Kahan summation temps for f

    vector<int> edge;
    bool onedge = false;

    for (it = circumtri.begin(); it != circumtri.end(); it++) {
        int t = *it;
        double vx = INDEX2(this->centers, t, 0);
        double vy = INDEX2(this->centers, t, 1);
        vector<double> c(6);
        //vector<bool> onedge(3);
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
                    onedge = true;
                    edge.push_back(INDEX3(this->nodes, t, j));
                    edge.push_back(INDEX3(this->nodes, t, k));
                }
            }
        }
        for (int i=0; i<3; i++) {
            int j = EDGE0(i);
            int k = EDGE1(i);
            int q = INDEX3(this->nodes, t, i);
            double ati = 0.0;

            if (!onedge || ((edge[0] != q) && edge[1] != q)) {
                ati = signed_area(vx, vy, 
                                  INDEX2(c, j, 0), INDEX2(c, j, 1),
                                  INDEX2(c, k, 0), INDEX2(c, k, 1));

                yA = ati - cA;
                tA = A + yA;
                cA = (tA - A) - yA;
                A = tA;

                yf = ati*z[q] - cf;
                tf = f + yf;
                cf = (tf - f) - yf;
                f = tf;
            }
        }
    }

    // If we're on an edge, then the scheme of adding up triangles as above
    // doesn't work so well. We'll take care of these two nodes here.
    if (onedge) {
        if (test) cout << "We're on an edge! "<< targetx<<" "<<targety << endl;
        set<int> T(circumtri.begin(), circumtri.end());
        vector<int> newedges0; // the two nodes that edge[0] still connect to
        vector<int> newedges1; // the two nodes that edge[1] still connect to
        set<int> alltri0; // all of the circumtriangle edge[0] participates in
        set<int> alltri1; // all of the circumtriangle edge[1] participates in
        for (it = circumtri.begin(); it != circumtri.end(); it++) {
            for (int i=0; i<3; i++) {
                int ti = INDEX3(this->neighbors, *it, i);
                int j = EDGE0(i);
                int k = EDGE1(i);
                int q0 = INDEX3(this->nodes, *it, j);
                int q1 = INDEX3(this->nodes, *it, k);

                if ((q0 == edge[0]) || (q1 == edge[0])) alltri0.insert(*it);
                if ((q0 == edge[1]) || (q1 == edge[1])) alltri1.insert(*it);

                if (!T.count(ti)) {
                    // neighbor is not in the set of circumtriangles
                    if (q0 == edge[0]) newedges0.push_back(q1);
                    if (q1 == edge[0]) newedges0.push_back(q0);
                    if (q0 == edge[1]) newedges1.push_back(q1);
                    if (q1 == edge[1]) newedges1.push_back(q0);
                }
            }
        }

        if (test) {
        cout << "  newedges0 ";
        copy(newedges0.begin(), newedges0.end(), ostream_iterator<int>(cout, " "));
        cout << endl;
        cout << "  newedges1 ";
        copy(newedges1.begin(), newedges1.end(), ostream_iterator<int>(cout, " "));
        cout << endl;
        cout << "  alltri0 ";
        copy(alltri0.begin(), alltri0.end(), ostream_iterator<int>(cout, " "));
        cout << endl;
        cout << "  alltri1 ";
        copy(alltri1.begin(), alltri1.end(), ostream_iterator<int>(cout, " "));
        cout << endl;
        }

        double cx, cy;
        circumcenter(this->x[edge[0]], this->y[edge[0]],
                     this->x[newedges0[0]], this->y[newedges0[0]],
                     targetx, targety,
                     cx, cy);
        if (test) cout << cx <<","<<cy<<endl;
        ConvexPolygon poly0(cx, cy);
        circumcenter(this->x[edge[0]], this->y[edge[0]],
                     this->x[newedges0[1]], this->y[newedges0[1]],
                     targetx, targety,
                     cx, cy);
        if (test) cout << cx <<","<<cy<<endl;
        poly0.push(cx, cy);

        circumcenter(this->x[edge[1]], this->y[edge[1]],
                     this->x[newedges1[0]], this->y[newedges1[0]],
                     targetx, targety,
                     cx, cy);
        //if (test) cout << cx <<","<<cy<<endl;
        ConvexPolygon poly1(cx, cy);
        circumcenter(this->x[edge[1]], this->y[edge[1]],
                     this->x[newedges1[1]], this->y[newedges1[1]],
                     targetx, targety,
                     cx, cy);
        //if (test) cout << cx <<","<<cy<<endl;
        poly1.push(cx, cy);

        set<int>::iterator sit;
        for (sit = alltri0.begin(); sit != alltri0.end(); sit++) {
            poly0.push(INDEX2(this->centers, *sit, 0),
                       INDEX2(this->centers, *sit, 1)); 
            if (test) cout <<*sit<<" "<<INDEX2(this->centers, *sit, 0)<<","<<
                       INDEX2(this->centers, *sit, 1) <<endl;
        }
        for (sit = alltri1.begin(); sit != alltri1.end(); sit++) {
            poly1.push(INDEX2(this->centers, *sit, 0),
                       INDEX2(this->centers, *sit, 1)); 
        }

        double a0 = poly0.area();
        double a1 = poly1.area();

        if (test) {
            cout << "Edge 0 "<<edge[0]<<" "<<x[edge[0]]<<","<<y[edge[0]]<<endl;
            cout << "Edge 1 "<<edge[1]<<" "<<x[edge[1]]<<","<<y[edge[1]]<<endl;
            cout << "  a0 = " << a0 <<endl;
            cout << "  a1 = " << a1 <<endl;

            vector<SeededPoint>::iterator spit;

            for (spit=poly0.points.begin(); spit!=poly0.points.end(); spit++) {
                cout << "  "<< spit->x <<","<<spit->y;
            }
            cout << "  " << poly0.x0<<","<<poly0.y0<<endl;

            for (spit=poly1.points.begin(); spit!=poly1.points.end(); spit++) {
                cout << "  "<<spit->x<<","<<spit->y;
            }
            cout << "  " << poly1.x0<<","<<poly1.y0<<endl;
        }

        f += a0*z[edge[0]];
        A += a0;
        f += a1*z[edge[1]];
        A += a1;

        // Anticlimactic, isn't it?
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

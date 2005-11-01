
#ifndef _DELAUNAY_UTILS_H
#define _DELAUNAY_UTILS_H

#define ONRIGHT(x0, y0, x1, y1, x, y) ((y0-y)*(x1-x) > (x0-x)*(y1-y))
#define EDGE0(node) ((node + 1) % 3)
#define EDGE1(node) ((node + 2) % 3)
#define INDEX2(arr,ix,jx) (arr[2*ix+jx])
#define INDEX3(arr,ix,jx) (arr[3*ix+jx])
#define INDEXN(arr,N,ix,jx) (arr[N*ix+jx])
#define SQ(a) ((a)*(a))

#define TOLERANCE_EPS (2e-16)
#define PERTURB_EPS (16*TOLERANCE_EPS)

extern int walking_triangles(int start, double targetx, double targety, 
    double *x, double *y, int *nodes, int *neighbors);
extern void getminmax(double *arr, int n, double& minimum, double& maximum);
extern bool circumcenter(double x0, double y0,
                         double x1, double y1,
                         double x2, double y2,
                         double& centerx, double& centery);
extern double signed_area(double x0, double y0,
                          double x1, double y1,
                          double x2, double y2);


#endif // _DELAUNAY_UTILS_H

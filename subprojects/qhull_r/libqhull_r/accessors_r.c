#include <stdlib.h>
#include "libqhull_r.h"

qhT *qh_alloc_qh(FILE *errfile) {
  qhT *qh = (qhT*) qh_malloc(sizeof(qhT));
  QHULL_LIB_CHECK
  if (qh) qh_zero(qh, errfile);
  return qh;
}

/* Free the qhT pointer.

   Note: qh_freeqhull and qh_memfreeshort should be called before calling qh_free_qh. */
void qh_free_qh(qhT *qh) {
  qh_free(qh);
}

#define GETTER(TYPE, FIELD) TYPE qh_get_##FIELD(const qhT *qh) { return qh->FIELD; }
#define SETTER(TYPE, FIELD) void qh_set_##FIELD(qhT *qh, TYPE _val_) { qh->FIELD = _val_; }

/* required to emulate the various FOR* macros */
GETTER(facetT*, facet_list)
GETTER(pointT*, first_point)
GETTER(int, hull_dim)
GETTER(int, num_facets)
GETTER(int, num_points)
GETTER(int, num_vertices)
GETTER(vertexT*, vertex_list)

/* suggested by @blegat based on Qhull.jl wrappers */
GETTER(realT, totarea)
GETTER(realT, totvol)
GETTER(boolT, hasAreaVolume)
SETTER(boolT, hasAreaVolume)
GETTER(boolT, hasTriangulation)
SETTER(boolT, hasTriangulation)

/* suggested by @JuhaHeiskala based on his DirectQHull.jl wrapper */
GETTER(int, num_good)
GETTER(setT*, del_vertices)
GETTER(int, input_dim)

/* other accessors, mimicking those provided by the scipy API: */
GETTER(boolT, DELAUNAY)
GETTER(boolT, SCALElast)
GETTER(boolT, KEEPcoplanar)
GETTER(boolT, MERGEexact)
GETTER(boolT, NOerrexit)
GETTER(boolT, PROJECTdelaunay)
GETTER(boolT, ATinfinity)
GETTER(boolT, UPPERdelaunay)
GETTER(int, normal_size)
GETTER(int, num_visible)
GETTER(int, center_size)
GETTER(const char *, qhull_command)
GETTER(facetT*, facet_tail)
GETTER(vertexT*, vertex_tail)
GETTER(unsigned int, facet_id)
GETTER(unsigned int, visit_id)
GETTER(unsigned int, vertex_visit)
GETTER(pointT*, input_points)
GETTER(coordT*, feasible_point)
GETTER(realT, last_low)
GETTER(realT, last_high)
GETTER(realT, last_newhigh)
GETTER(realT, max_outside)
GETTER(realT, MINoutside)
GETTER(realT, DISTround)
GETTER(setT*, other_points)

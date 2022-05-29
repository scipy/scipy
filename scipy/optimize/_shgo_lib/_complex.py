"""
Base classes for low memory simplicial complex structures.
"""
# Std. Library
import copy
import collections
import logging
import os
import itertools
import json
import decimal
from functools import lru_cache  # For Python 3 only
from abc import ABC, abstractmethod
# Required modules:
import numpy

# Module specific imports
from ._vertex import VertexCacheField

import matplotlib
from matplotlib import pyplot
from matplotlib.patches import FancyArrowPatch
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d import axes3d, Axes3D, proj3d

matplotlib_available = True



# Main complex class:
class Complex:
    def __init__(self, dim, domain=None, sfield=None, sfield_args=(),
                 symmetry=None, constraints=None, workers=1):
        """
        A base class for a simplicial complex described as a cache of vertices
        together with their connections.

        Important methods:
            Domain triangulation:
                    Complex.triangulate, Complex.split_generation
            Triangulating arbitrary points (must be traingulable,
                may exist outside domain):
                    Complex.triangulate(sample_set)
            Converting another simplicial complex structure data type to the
                structure used in Complex (ex. OBJ wavefront)
                    Complex.convert(datatype, data)

        Important objects:
            HC.V: The cache of vertices and their connection
            HC.H: Storage structure of all vertex groups

        :param dim: int, Spatial dimensionality of the complex R^dim
        :param domain: list of tuples, optional
                The bounds [x_l, x_u]^dim of the hyperrectangle space
                ex. The default domain is the hyperrectangle [0, 1]^dim
                Note: The domain must be convex, non-convex spaces can be cut
                      away from this domain using the non-linear
                      g_cons functions to define any arbitrary domain
                      (these domains may also be disconnected from each other)
        :param sfield: A scalar function defined in the associated domain
                           f: R^dim --> R
        :param sfield_args: tuple, Additional arguments to be passed to sfield
        :param vfield: A scalar function defined in the associated domain
                           f: R^dim --> R^m
                       (for example a gradient function of the scalar field)
        :param vfield_args: tuple, Additional arguments to be passed to sfield
        :param symmetry: If all the variables in the field are symmetric this
                option will reduce complexity of the triangulation by O(n!)

        constraints : dict or sequence of dict, optional
            Constraints definition.
            Function(s) ``R**n`` in the form::

                g(x) <= 0 applied as g : R^n -> R^m
                h(x) == 0 applied as h : R^n -> R^p

            Each constraint is defined in a dictionary with fields:

                type : str
                    Constraint type: 'eq' for equality, 'ineq' for inequality.
                fun : callable
                    The function defining the constraint.
                jac : callable, optional
                    The Jacobian of `fun` (only for SLSQP).
                args : sequence, optional
                    Extra arguments to be passed to the function and Jacobian.

            Equality constraint means that the constraint function result is to
            be zero whereas inequality means that it is to be non-negative.constraints : dict or sequence of dict, optional
            Constraints definition.
            Function(s) ``R**n`` in the form::

                g(x) <= 0 applied as g : R^n -> R^m
                h(x) == 0 applied as h : R^n -> R^p

            Each constraint is defined in a dictionary with fields:

                type : str
                    Constraint type: 'eq' for equality, 'ineq' for inequality.
                fun : callable
                    The function defining the constraint.
                jac : callable, optional
                    The Jacobian of `fun` (unused).
                args : sequence, optional
                    Extra arguments to be passed to the function and Jacobian.

            Equality constraint means that the constraint function result is to
            be zero whereas inequality means that it is to be non-negative.

        workers : int  optional
            Uses `multiprocessing.Pool <multiprocessing>`) to compute the field
             functions in parrallel.
        """
        self.dim = dim

        # Domains
        self.domain = domain
        if domain is None:
            self.bounds = [(float(0), float(1.0)), ] * dim
        else:
            self.bounds = domain
        self.symmetry = symmetry
        #      here in init to avoid if checks

        # Field functions
        self.sfield = sfield
        self.sfield_args = sfield_args

        # Process constraints
        # Constraints
        # Process constraint dict sequence:
        if constraints is not None:
            self.min_cons = constraints
            self.g_cons = []
            self.g_args = []
            if (type(constraints) is not tuple) and (type(constraints)
                                                     is not list):
                constraints = (constraints,)

            for cons in constraints:
                if cons['type'] == 'ineq':
                    self.g_cons.append(cons['fun'])
                    try:
                        self.g_args.append(cons['args'])
                    except KeyError:
                        self.g_args.append(())
            self.g_cons = tuple(self.g_cons)
            self.g_args = tuple(self.g_args)
        else:
            self.g_cons = None
            self.g_args = None

        # Homology properties
        self.gen = 0
        self.perm_cycle = 0

        # Every cell is stored in a list of its generation,
        # ex. the initial cell is stored in self.H[0]
        # 1st get new cells are stored in self.H[1] etc.
        # When a cell is sub-generated it is removed from this list

        self.H = []  # Storage structure of vertex groups

        # Cache of all vertices
        if (sfield is not None) or (self.g_cons is not None):
            # Initiate a vertex cache and an associated field cache, note that
            # the field case is always initiated inside the vertex cache if an
            # associated field scalar field is defined:
            if sfield is not None:
                self.V = VertexCacheField(field=sfield, field_args=sfield_args,
                                          g_cons=self.g_cons,
                                          g_cons_args=self.g_args,
                                          workers=workers)
            elif self.g_cons is not None:
                self.V = VertexCacheField(field=sfield, field_args=sfield_args,
                                          g_cons=self.g_cons,
                                          g_cons_args=self.g_args,
                                          workers=workers)
        else:
            self.V = VertexCacheIndex()

        self.V_non_symm = []  # List of non-symmetric vertices

    def __call__(self):
        return self.H

    # %% Triangulation methods
    def cyclic_product(self, bounds, origin, supremum, centroid=True):
        """Generate initial triangulation using cyclic product"""

        vot = tuple(origin)
        vut = tuple(supremum)  # Hyperrectangle supremum
        self.V[vot]
        vo = self.V[vot]
        yield vo.x
        self.V[vut].connect(self.V[vot])
        yield vut
        # Cyclic group approach with second x_l --- x_u operation.

        # These containers store the "lower" and "upper" vertices
        # corresponding to the origin or supremum of every C2 group.
        # It has the structure of `dim` times embedded lists each containing
        # these vertices as the entire complex grows. Bounds[0] has to be done
        # outside the loops before we have symmetric containers.
        #NOTE: This means that bounds[0][1] must always exist
        C0x = [[self.V[vot]]]
        a_vo = copy.copy(list(origin))
        a_vo[0] = vut[0]  # Update aN Origin
        a_vo = self.V[tuple(a_vo)]
        #self.V[vot].connect(self.V[tuple(a_vo)])
        self.V[vot].connect(a_vo)
        yield a_vo.x
        C1x = [[a_vo]]
        #C1x = [[self.V[tuple(a_vo)]]]
        ab_C = []  # Container for a + b operations

        # Loop over remaining bounds
        for i, x in enumerate(bounds[1:]):
            # Update lower and upper containers
            C0x.append([])
            C1x.append([])
            # try to access a second bound (if not, C1 is symmetric)
            try:
                # Early try so that we don't have to copy the cache before
                # moving on to next C1/C2: Try to add the operation of a new
                # C2 product by accessing the upper bound
                x[1]
                # Copy lists for iteration
                cC0x = [x[:] for x in C0x[:i + 1]]
                cC1x = [x[:] for x in C1x[:i + 1]]
                for j, (VL, VU) in enumerate(zip(cC0x, cC1x)):
                    for k, (vl, vu) in enumerate(zip(VL, VU)):
                        # Build aN vertices for each lower-upper pair in N:
                        a_vl = list(vl.x)
                        a_vu = list(vu.x)
                        a_vl[i + 1] = vut[i + 1]
                        a_vu[i + 1] = vut[i + 1]
                        a_vl = self.V[tuple(a_vl)]

                        # Connect vertices in N to corresponding vertices
                        # in aN:
                        vl.connect(a_vl)

                        yield a_vl.x

                        a_vu = self.V[tuple(a_vu)]
                        # Connect vertices in N to corresponding vertices
                        # in aN:
                        vu.connect(a_vu)

                        # Connect new vertex pair in aN:
                        a_vl.connect(a_vu)

                        # Connect lower pair to upper (triangulation
                        # operation of a + b (two arbitrary operations):
                        vl.connect(a_vu)
                        ab_C.append((vl, a_vu))

                        # Update the containers
                        C0x[i + 1].append(vl)
                        C0x[i + 1].append(vu)
                        C1x[i + 1].append(a_vl)
                        C1x[i + 1].append(a_vu)

                        # Update old containers
                        C0x[j].append(a_vl)
                        C1x[j].append(a_vu)

                        # Yield new points
                        yield a_vu.x

                # Try to connect aN lower source of previous a + b
                # operation with a aN vertex
                ab_Cc = copy.copy(ab_C)

                #TODO: SHOULD THIS BE MOVED OUTSIDE THE try ?
                for vp in ab_Cc:
                    b_v = list(vp[0].x)
                    ab_v = list(vp[1].x)
                    b_v[i + 1] = vut[i + 1]
                    ab_v[i + 1] = vut[i + 1]
                    b_v = self.V[tuple(b_v)]  # b + vl
                    ab_v = self.V[tuple(ab_v)]  # b + a_vl
                    # Note o---o is already connected
                    vp[0].connect(ab_v)  # o-s
                    b_v.connect(ab_v)  # s-s

                    # Add new list of cross pairs
                    ab_C.append((vp[0], ab_v))
                    ab_C.append((b_v, ab_v))

            except IndexError:
                cC0x = C0x[i]
                cC1x = C1x[i]
                VL, VU = cC0x, cC1x
                for k, (vl, vu) in enumerate(zip(VL, VU)):
                    # Build aN vertices for each lower-upper pair in N:
                    a_vu = list(vu.x)
                    a_vu[i + 1] = vut[i + 1]
                    # Connect vertices in N to corresponding vertices
                    # in aN:
                    a_vu = self.V[tuple(a_vu)]
                    # Connect vertices in N to corresponding vertices
                    # in aN:
                    vu.connect(a_vu)
                    # Connect new vertex pair in aN:
                    # a_vl.connect(a_vu)
                    # Connect lower pair to upper (triangulation
                    # operation of a + b (two arbitrary operations):
                    vl.connect(a_vu)
                    ab_C.append((vl, a_vu))
                    C0x[i + 1].append(vu)
                    C1x[i + 1].append(a_vu)
                    # Yield new points
                    a_vu.connect(self.V[vut])
                    yield a_vu.x
                    ab_Cc = copy.copy(ab_C)
                    for vp in ab_Cc:
                        if vp[1].x[i] == vut[i]:
                            ab_v = list(vp[1].x)
                            ab_v[i + 1] = vut[i + 1]
                            ab_v = self.V[tuple(ab_v)]  # b + a_vl
                            # Note o---o is already connected
                            vp[0].connect(ab_v)  # o-s

                            # Add new list of cross pairs
                            ab_C.append((vp[0], ab_v))


        # Clean class trash
        try:
            del C0x
            del cC0x
            del C1x
            del cC1x
            del ab_C
            del ab_Cc
        except UnboundLocalError:
            pass

        # Extra yield to ensure that the triangulation is completed
        if centroid:
            vo = self.V[vot]
            vs = self.V[vut]
            # Disconnect the origin and supremum
            vo.disconnect(vs)
            # Build centroid
            vc = self.split_edge(vot, vut)
            for v in vo.nn:
                v.connect(vc)
            yield vc.x
            return vc.x
        else:
            yield vut
            return vut


    def triangulate(self, n=None, symmetry=None, centroid=True, printout=False):
        """
        Triangulate the initial domain, if n is not None then a limited number
        of points will be generated

        :param n:
        :param symmetry:

            Ex. Dictionary/hashtable
            f(x) = (x_1 + x_2 + x_3) + (x_4)**2 + (x_5)**2 + (x_6)**2

            symmetry = symmetry[0]: 0,  # Variable 1
                       symmetry[1]: 0,  # symmetric to variable 1
                       symmetry[2]: 0,  # symmetric to variable 1
                       symmetry[3]: 3,  # Variable 4
                       symmetry[4]: 3,  # symmetric to variable 4
                       symmetry[5]: 3,  # symmetric to variable 4
                        }

        :param printout:
        :return:

        NOTES:
        ------
        Rather than using the combinatorial algorithm to connect vertices we
        make the following observation:

        The bound pairs are similar a C2 cyclic group and the structure is
        formed using the cartesian product:

        H = C2 x C2 x C2 ... x C2 (dim times)

        So construct any normal subgroup N and consider H/N first, we connect
        all vertices within N (ex. N is C2 (the first dimension), then we move
        to a left coset aN (an operation moving around the defined H/N group by
        for example moving from the lower bound in C2 (dimension 2) to the
        higher bound in C2. During this operation connection all the vertices.
        Now repeat the N connections. Note that these elements can be connected
        in parrallel.
        """
        # Inherit class arguments
        if symmetry is None:
            symmetry = self.symmetry
        # Build origin and supremum vectors
        origin = [i[0] for i in self.bounds]
        self.origin = origin
        supremum = [i[1] for i in self.bounds]

        self.supremum = supremum

        if symmetry is None:
            cbounds = self.bounds
        else:
            cbounds = copy.copy(self.bounds)
            for i, j in enumerate(symmetry):
                if i is not j:
                    # pop second entry on second symmetry vars
                    cbounds[i] = [self.bounds[symmetry[i]][0]]
                    # Sole (first) entry is the sup value and there is no origin
                    cbounds[i] = [self.bounds[symmetry[i]][1]]
                    if self.bounds[symmetry[i]] is not self.bounds[symmetry[j]]:
                        logging.warning(f"Variable {i} was specified as "
                                        f"symmetetric to variable {j}, however,"
                                        f"the bounds {i} ="
                                        f" {self.bounds[symmetry[i]]} and {j} ="
                                        f" {self.bounds[symmetry[j]]} do not "
                                        f"match, the mismatch was ignored in "
                                        f"the initial triangulation.")
                        cbounds[i] = self.bounds[symmetry[j]]

        if n is None:
            # Build generator
            self.cp = self.cyclic_product(cbounds, origin, supremum, centroid)
            for i in self.cp:
                i

            try:
                self.triangulated_vectors.append((tuple(self.origin),
                                                  tuple(self.supremum)))
            except (AttributeError, KeyError):
                self.triangulated_vectors = [(tuple(self.origin),
                                              tuple(self.supremum))]

        else:
            # Check if generator already exists
            try:
                self.cp
            except (AttributeError, KeyError):
                self.cp = self.cyclic_product(cbounds, origin, supremum,
                                              centroid)

            try:
                while len(self.V.cache) < n:
                    next(self.cp)
            except StopIteration:
                try:
                    self.triangulated_vectors.append((tuple(self.origin),
                                                      tuple(self.supremum)))
                except (AttributeError, KeyError):
                    self.triangulated_vectors = [(tuple(self.origin),
                                                  tuple(self.supremum))]

        if printout:
            # for v in self.C0():
            #   v.print_out()
            for v in self.V.cache:
                self.V[v].print_out()

        return

    def refine(self, n=1):
        if n is None:
            try:
                self.triangulated_vectors
                self.refine_all()
                return
            except AttributeError as ae:
                if str(ae) == "'Complex' object has no attribute " \
                              "'triangulated_vectors'":
                    self.triangulate(symmetry=self.symmetry)
                    return
                else:
                    raise

        nt = len(self.V.cache) + n  # Target number of total vertices
        # In the outer while loop we iterate until we have added an extra `n`
        # vertices to the complex:
        while len(self.V.cache) < nt:  # while loop 1
            try:  # try 1
                # Try to access triangulated_vectors, this should only be
                # defined if an initial triangulation has already been performed
                self.triangulated_vectors
                # Try a usual iteration of the current generator, if it
                # does not exist or is exhausted then produce a new generator
                try:  # try 2
                    next(self.rls)
                except (AttributeError, StopIteration, KeyError):
                    vp = self.triangulated_vectors[0]
                    self.rls = self.refine_local_space(*vp, bounds=self.bounds)
                    next(self.rls)

            except (AttributeError, KeyError):
                # If an initial triangulation has not been completed, then
                # we start/continue the initial triangulation targeting `nt`
                # vertices, if nt is greater than the initial number of vertices
                # then the refine routine will move back to try 1.
                self.triangulate(nt, self.symmetry)
        return

    def refine_all(self, centroids=True):
        """
        Refine the entire domain of the current complex
        :return:
        """
        try:
            self.triangulated_vectors
            tvs = copy.copy(self.triangulated_vectors)
            for i, vp in enumerate(tvs):
                self.rls = self.refine_local_space(*vp, bounds=self.bounds)
                for i in self.rls:
                    i
        except AttributeError as ae:
            if str(ae) == "'Complex' object has no attribute " \
                          "'triangulated_vectors'":
                self.triangulate(symmetry=self.symmetry, centroid=centroids)
            else:
                raise

        # This adds a centroid to every new sub-domain generated and defined
        # by self.triangulated_vectors, in addition the vertices ! to complete
        # the triangulation
        return

    def refine_local_space(self, origin, supremum, bounds, centroid=1):
        # Copy for later removal
        origin_c = copy.copy(origin)
        supremum_c = copy.copy(supremum)

        # Change the vector orientation so that it is only increasing
        s_ov = list(origin)
        s_origin = list(origin)
        s_sv = list(supremum)
        s_supremum = list(supremum)
        for i, vi in enumerate(s_origin):
            if s_ov[i] > s_sv[i]:
                s_origin[i] = s_sv[i]
                s_supremum[i] = s_ov[i]

        vot = tuple(s_origin)
        vut = tuple(s_supremum)  # Hyperrectangle supremum

        vo = self.V[vot]  # initiate if doesn't exist yet
        vs = self.V[vut]
        # Start by finding the old centroid of the new space:
        vco = self.split_edge(vo.x, vs.x)  # Split in case not centroid arg

        # Find set of extreme vertices in current local space
        sup_set = copy.copy(vco.nn)
        # Cyclic group approach with second x_l --- x_u operation.

        # These containers store the "lower" and "upper" vertices
        # corresponding to the origin or supremum of every C2 group.
        # It has the structure of `dim` times embedded lists each containing
        # these vertices as the entire complex grows. Bounds[0] has to be done
        # outside the loops before we have symmetric containers.
        #NOTE: This means that bounds[0][1] must always exist

        a_vl = copy.copy(list(vot))
        a_vl[0] = vut[0]  # Update aN Origin
        if tuple(a_vl) not in self.V.cache:
            vo = self.V[vot]  # initiate if doesn't exist yet
            vs = self.V[vut]
            # Start by finding the old centroid of the new space:
            vco = self.split_edge(vo.x, vs.x)  # Split in case not centroid arg

            # Find set of extreme vertices in current local space
            sup_set = copy.copy(vco.nn)
            a_vl = copy.copy(list(vot))
            a_vl[0] = vut[0]  # Update aN Origin
            a_vl = self.V[tuple(a_vl)]
        else:
            a_vl = self.V[tuple(a_vl)]

        c_v = self.split_edge(vo.x, a_vl.x)
        c_v.connect(vco)
        yield c_v.x
        Cox = [[vo]]
        Ccx = [[c_v]]
        Cux = [[a_vl]]
        ab_C = []  # Container for a + b operations
        s_ab_C = []  # Container for symmetric a + b operations

        # Loop over remaining bounds
        for i, x in enumerate(bounds[1:]):
            # Update lower and upper containers
            Cox.append([])
            Ccx.append([])
            Cux.append([])
            # try to access a second bound (if not, C1 is symmetric)
            try:
                t_a_vl = list(vot)
                t_a_vl[i + 1] = vut[i + 1]

                # New: lists are used anyway, so copy all
                # %%
                # Copy lists for iteration
                cCox = [x[:] for x in Cox[:i + 1]]
                cCcx = [x[:] for x in Ccx[:i + 1]]
                cCux = [x[:] for x in Cux[:i + 1]]
                # Try to connect aN lower source of previous a + b
                # operation with a aN vertex
                ab_Cc = copy.copy(ab_C)  # NOTE: We append ab_C in the
                # (VL, VC, VU) for-loop, but we use the copy of the list in the
                # ab_Cc for-loop.
                s_ab_Cc = copy.copy(s_ab_C)

                #%%
                # Early try so that we don't have to copy the cache before
                # moving on to next C1/C2: Try to add the operation of a new
                # C2 product by accessing the upper bound
                if tuple(t_a_vl) not in self.V.cache:
                    raise IndexError  # Raise error to continue symmetric refine
                t_a_vu = list(vut)
                t_a_vu[i + 1] = vut[i + 1]
                if tuple(t_a_vu) not in self.V.cache:
                    raise IndexError  # Raise error to continue symmetric refine

                for vectors in s_ab_Cc:
                    # s_ab_C.append([c_vc, vl, vu, a_vu])
                    bc_vc = list(vectors[0].x)
                    b_vl = list(vectors[1].x)
                    b_vu = list(vectors[2].x)
                    ba_vu = list(vectors[3].x)

                    bc_vc[i + 1] = vut[i + 1]
                    b_vl[i + 1] = vut[i + 1]
                    b_vu[i + 1] = vut[i + 1]
                    ba_vu[i + 1] = vut[i + 1]

                    bc_vc = self.V[tuple(bc_vc)]
                    bc_vc.connect(vco)  # NOTE: Unneeded?
                    yield bc_vc

                    # Split to centre, call this centre group "d = 0.5*a"
                    d_bc_vc = self.split_edge(vectors[0].x, bc_vc.x)
                    d_bc_vc.connect(bc_vc)
                    d_bc_vc.connect(vectors[1])  # Connect all to centroid
                    d_bc_vc.connect(vectors[2])  # Connect all to centroid
                    d_bc_vc.connect(vectors[3])  # Connect all to centroid
                    yield d_bc_vc.x
                    b_vl = self.V[tuple(b_vl)]
                    bc_vc.connect(b_vl)  # Connect aN cross pairs
                    d_bc_vc.connect(b_vl)  # Connect all to centroid

                    yield b_vl
                    b_vu = self.V[tuple(b_vu)]
                    bc_vc.connect(b_vu)  # Connect aN cross pairs
                    d_bc_vc.connect(b_vu)  # Connect all to centroid

                    b_vl_c = self.split_edge(b_vu.x, b_vl.x)
                    bc_vc.connect(b_vl_c)

                    yield b_vu
                    ba_vu = self.V[tuple(ba_vu)]
                    bc_vc.connect(ba_vu)  # Connect aN cross pairs
                    d_bc_vc.connect(ba_vu)  # Connect all to centroid

                    # Split the a + b edge of the initial triangulation:
                    os_v = self.split_edge(vectors[1].x, ba_vu.x)  # o-s
                    ss_v = self.split_edge(b_vl.x, ba_vu.x)  # s-s
                    b_vu_c = self.split_edge(b_vu.x, ba_vu.x)
                    bc_vc.connect(b_vu_c)
                    yield os_v.x  # often equal to vco, but not always
                    yield ss_v.x  # often equal to bc_vu, but not always
                    yield ba_vu
                    # Split remaining to centre, call this centre group "d = 0.5*a"
                    d_bc_vc = self.split_edge(vectors[0].x, bc_vc.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    yield d_bc_vc.x
                    d_b_vl = self.split_edge(vectors[1].x, b_vl.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    d_bc_vc.connect(d_b_vl)  # Connect dN cross pairs
                    yield d_b_vl.x
                    d_b_vu = self.split_edge(vectors[2].x, b_vu.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    d_bc_vc.connect(d_b_vu)  # Connect dN cross pairs
                    yield d_b_vu.x
                    d_ba_vu = self.split_edge(vectors[3].x, ba_vu.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    d_bc_vc.connect(d_ba_vu)  # Connect dN cross pairs
                    yield d_ba_vu

                    # comb = [c_vc, vl, vu, a_vl, a_vu,
                    #       bc_vc, b_vl, b_vu, ba_vl, ba_vu]
                    comb = [vl, vu, a_vu,
                            b_vl, b_vu, ba_vu]
                    comb_iter = itertools.combinations(comb, 2)
                    for vecs in comb_iter:
                        self.split_edge(vecs[0].x, vecs[1].x)
                    # Add new list of cross pairs
                    ab_C.append((d_bc_vc, vectors[1], b_vl, a_vu, ba_vu))
                    ab_C.append((d_bc_vc, vl, b_vl, a_vu, ba_vu))  # = prev

                for vectors in ab_Cc:
                    bc_vc = list(vectors[0].x)
                    b_vl = list(vectors[1].x)
                    b_vu = list(vectors[2].x)
                    ba_vl = list(vectors[3].x)
                    ba_vu = list(vectors[4].x)
                    bc_vc[i + 1] = vut[i + 1]
                    b_vl[i + 1] = vut[i + 1]
                    b_vu[i + 1] = vut[i + 1]
                    ba_vl[i + 1] = vut[i + 1]
                    ba_vu[i + 1] = vut[i + 1]
                    bc_vc = self.V[tuple(bc_vc)]
                    bc_vc.connect(vco)  # NOTE: Unneeded?
                    yield bc_vc

                    # Split to centre, call this centre group "d = 0.5*a"
                    d_bc_vc = self.split_edge(vectors[0].x, bc_vc.x)
                    d_bc_vc.connect(bc_vc)
                    d_bc_vc.connect(vectors[1])  # Connect all to centroid
                    d_bc_vc.connect(vectors[2])  # Connect all to centroid
                    d_bc_vc.connect(vectors[3])  # Connect all to centroid
                    d_bc_vc.connect(vectors[4])  # Connect all to centroid
                    yield d_bc_vc.x
                    b_vl = self.V[tuple(b_vl)]
                    bc_vc.connect(b_vl)  # Connect aN cross pairs
                    d_bc_vc.connect(b_vl)  # Connect all to centroid
                    yield b_vl
                    b_vu = self.V[tuple(b_vu)]
                    bc_vc.connect(b_vu)  # Connect aN cross pairs
                    d_bc_vc.connect(b_vu)  # Connect all to centroid
                    yield b_vu
                    ba_vl = self.V[tuple(ba_vl)]
                    bc_vc.connect(ba_vl)  # Connect aN cross pairs
                    d_bc_vc.connect(ba_vl)  # Connect all to centroid
                    self.split_edge(b_vu.x, ba_vl.x)
                    yield ba_vl
                    ba_vu = self.V[tuple(ba_vu)]
                    bc_vc.connect(ba_vu)  # Connect aN cross pairs
                    d_bc_vc.connect(ba_vu)  # Connect all to centroid
                    # Split the a + b edge of the initial triangulation:
                    os_v = self.split_edge(vectors[1].x, ba_vu.x)  # o-s
                    ss_v = self.split_edge(b_vl.x, ba_vu.x)  # s-s
                    yield os_v.x  # often equal to vco, but not always
                    yield ss_v.x  # often equal to bc_vu, but not always
                    yield ba_vu
                    # Split remaining to centre, call this centre group "d = 0.5*a"
                    d_bc_vc = self.split_edge(vectors[0].x, bc_vc.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    yield d_bc_vc.x
                    d_b_vl = self.split_edge(vectors[1].x, b_vl.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    d_bc_vc.connect(d_b_vl)  # Connect dN cross pairs
                    yield d_b_vl.x
                    d_b_vu = self.split_edge(vectors[2].x, b_vu.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    d_bc_vc.connect(d_b_vu)  # Connect dN cross pairs
                    yield d_b_vu.x
                    d_ba_vl = self.split_edge(vectors[3].x, ba_vl.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    d_bc_vc.connect(d_ba_vl)  # Connect dN cross pairs
                    yield d_ba_vl
                    d_ba_vu = self.split_edge(vectors[4].x, ba_vu.x)
                    d_bc_vc.connect(vco)  # NOTE: Unneeded?
                    d_bc_vc.connect(d_ba_vu)  # Connect dN cross pairs
                    yield d_ba_vu
                    c_vc, vl, vu, a_vl, a_vu = vectors

                    comb = [vl, vu, a_vl, a_vu,
                            b_vl, b_vu, ba_vl, ba_vu]
                    comb_iter = itertools.combinations(comb, 2)
                    for vecs in comb_iter:
                        self.split_edge(vecs[0].x, vecs[1].x)

                    # Add new list of cross pairs
                    ab_C.append((bc_vc, b_vl, b_vu, ba_vl, ba_vu))
                    ab_C.append((d_bc_vc, d_b_vl, d_b_vu, d_ba_vl, d_ba_vu))
                    ab_C.append((d_bc_vc, vectors[1], b_vl, a_vu, ba_vu))
                    ab_C.append((d_bc_vc, vu, b_vu, a_vl, ba_vl))

                for j, (VL, VC, VU) in enumerate(zip(cCox, cCcx, cCux)):
                    for k, (vl, vc, vu) in enumerate(zip(VL, VC, VU)):
                        # Build aN vertices for each lower-upper C3 group in N:
                        a_vl = list(vl.x)
                        a_vu = list(vu.x)
                        a_vl[i + 1] = vut[i + 1]
                        a_vu[i + 1] = vut[i + 1]
                        a_vl = self.V[tuple(a_vl)]
                        a_vu = self.V[tuple(a_vu)]
                        # Note, build (a + vc) later for consistent yields
                        # Split the a + b edge of the initial triangulation:
                        c_vc = self.split_edge(vl.x, a_vu.x)
                        self.split_edge(vl.x, vu.x)  # Equal to vc
                        # Build cN vertices for each lower-upper C3 group in N:
                        c_vc.connect(vco)
                        c_vc.connect(vc)
                        c_vc.connect(vl)  # Connect c + ac operations
                        c_vc.connect(vu)  # Connect c + ac operations
                        c_vc.connect(a_vl)  # Connect c + ac operations
                        c_vc.connect(a_vu)  # Connect c + ac operations
                        yield(c_vc.x)
                        c_vl = self.split_edge(vl.x, a_vl.x)
                        c_vl.connect(vco)
                        c_vc.connect(c_vl)  # Connect cN group vertices
                        yield(c_vl.x)
                        c_vu = self.split_edge(vu.x, a_vu.x) # yield at end of loop
                        c_vu.connect(vco)
                        # Connect remaining cN group vertices
                        c_vc.connect(c_vu)  # Connect cN group vertices
                        yield (c_vu.x)

                        a_vc = self.split_edge(a_vl.x, a_vu.x)  # is (a + vc) ?
                        a_vc.connect(vco)
                        a_vc.connect(c_vc)

                        # Storage for connecting c + ac operations:
                        ab_C.append((c_vc, vl, vu, a_vl, a_vu))

                        # Update the containers
                        Cox[i + 1].append(vl)
                        Cox[i + 1].append(vc)
                        Cox[i + 1].append(vu)
                        Ccx[i + 1].append(c_vl)
                        Ccx[i + 1].append(c_vc)
                        Ccx[i + 1].append(c_vu)
                        Cux[i + 1].append(a_vl)
                        Cux[i + 1].append(a_vc)
                        Cux[i + 1].append(a_vu)

                        # Update old containers
                        Cox[j].append(c_vl)  # !
                        Cox[j].append(a_vl)
                        Ccx[j].append(c_vc)  # !
                        Ccx[j].append(a_vc)  # !
                        Cux[j].append(c_vu)  # !
                        Cux[j].append(a_vu)

                        # Yield new points
                        yield(a_vc.x)

            except IndexError:
                for vectors in ab_Cc:
                    ba_vl = list(vectors[3].x)
                    ba_vu = list(vectors[4].x)
                    ba_vl[i + 1] = vut[i + 1]
                    ba_vu[i + 1] = vut[i + 1]
                    ba_vu = self.V[tuple(ba_vu)]
                    yield ba_vu
                    d_bc_vc = self.split_edge(vectors[1].x, ba_vu.x)  # o-s
                    yield ba_vu
                    d_bc_vc.connect(vectors[1])  # Connect all to centroid
                    d_bc_vc.connect(vectors[2])  # Connect all to centroid
                    d_bc_vc.connect(vectors[3])  # Connect all to centroid
                    d_bc_vc.connect(vectors[4])  # Connect all to centroid
                    yield d_bc_vc.x
                    ba_vl = self.V[tuple(ba_vl)]
                    yield ba_vl
                    d_ba_vl = self.split_edge(vectors[3].x, ba_vl.x)
                    d_ba_vu = self.split_edge(vectors[4].x, ba_vu.x)
                    d_ba_vc = self.split_edge(d_ba_vl.x, d_ba_vu.x)
                    yield d_ba_vl
                    yield d_ba_vu
                    yield d_ba_vc
                    c_vc, vl, vu, a_vl, a_vu = vectors
                    comb = [vl, vu, a_vl, a_vu,
                            ba_vl,
                            ba_vu]
                    comb_iter = itertools.combinations(comb, 2)
                    for vecs in comb_iter:
                        self.split_edge(vecs[0].x, vecs[1].x)

                # Copy lists for iteration
                cCox = Cox[i]
                cCcx = Ccx[i]
                cCux = Cux[i]
                VL, VC, VU = cCox, cCcx, cCux
                for k, (vl, vc, vu) in enumerate(zip(VL, VC, VU)):
                    # Build aN vertices for each lower-upper pair in N:
                    a_vu = list(vu.x)
                    a_vu[i + 1] = vut[i + 1]

                    # Connect vertices in N to corresponding vertices
                    # in aN:
                    a_vu = self.V[tuple(a_vu)]
                    yield a_vl.x
                    # Split the a + b edge of the initial triangulation:
                    c_vc = self.split_edge(vl.x, a_vu.x)
                    self.split_edge(vl.x, vu.x)  # Equal to vc
                    c_vc.connect(vco)
                    c_vc.connect(vc)
                    c_vc.connect(vl)  # Connect c + ac operations
                    c_vc.connect(vu)  # Connect c + ac operations
                    c_vc.connect(a_vu)  # Connect c + ac operations
                    yield (c_vc.x)
                    c_vu = self.split_edge(vu.x,
                                           a_vu.x)  # yield at end of loop
                    c_vu.connect(vco)
                    # Connect remaining cN group vertices
                    c_vc.connect(c_vu)  # Connect cN group vertices
                    yield (c_vu.x)

                    # Update the containers
                    Cox[i + 1].append(vu)
                    Ccx[i + 1].append(c_vu)
                    Cux[i + 1].append(a_vu)

                    # Update old containers
                    s_ab_C.append([c_vc, vl, vu, a_vu])

                    yield a_vu.x

        # Clean class trash
        try:
            del C0x
            del cC0x
            del C1x
            del cC1x
            del ab_C
            del ab_Cc
        except UnboundLocalError:
            pass

        try:
            self.triangulated_vectors.remove((tuple(origin_c),
                                              tuple(supremum_c)))
        except ValueError:
            # Turn this into a logging warning?
            pass
        # Add newly triangulated vectors:
        for vs in sup_set:
            self.triangulated_vectors.append((tuple(vco.x), tuple(vs.x)))

        # Extra yield to ensure that the triangulation is completed
        if centroid:
            vcn_set = set()
            c_nn_lists = []
            for vs in sup_set:
                # Build centroid
                c_nn = self.vpool(vco.x, vs.x)
                try:
                    c_nn.remove(vcn_set)
                except KeyError:
                    pass
                c_nn_lists.append(c_nn)

            for c_nn in c_nn_lists:
                try:
                    c_nn.remove(vcn_set)
                except KeyError:
                    pass

            for vs, c_nn in zip(sup_set, c_nn_lists):
                # Build centroid
                vcn = self.split_edge(vco.x, vs.x)
                vcn_set.add(vcn)
                try:  # Shouldn't be needed?
                    c_nn.remove(vcn_set)
                except KeyError:
                    pass
                for vnn in c_nn:
                    vcn.connect(vnn)
                yield vcn.x
        else:
            pass

        yield vut
        return

    def refine_star(self, v):
        """
        Refine the star domain of a vertex v
        :param v: Vertex object
        """
        # Copy lists before iteration
        vnn = copy.copy(v.nn)
        v1nn = []
        d_v0v1_set = set()
        for v1 in vnn:
            v1nn.append(copy.copy(v1.nn))

        for v1, v1nn in zip(vnn, v1nn):
            vnnu = v1nn.intersection(vnn)

            d_v0v1 = self.split_edge(v.x, v1.x)
            for o_d_v0v1 in d_v0v1_set:
                d_v0v1.connect(o_d_v0v1)
            d_v0v1_set.add(d_v0v1)
            for v2 in vnnu:
                d_v1v2 = self.split_edge(v1.x, v2.x)
                d_v0v1.connect(d_v1v2)
        return

    @lru_cache(maxsize=None)
    def split_edge(self, v1, v2):
        v1 = self.V[v1]
        v2 = self.V[v2]
        # Destroy original edge, if it exists:
        v1.disconnect(v2)
        # Compute vertex on centre of edge:
        try:
            vct = (v2.x_a - v1.x_a) / 2.0 + v1.x_a
        except TypeError:  # Allow for decimal operations
            vct = (v2.x_a - v1.x_a) / decimal.Decimal(2.0) + v1.x_a

        vc = self.V[tuple(vct)]
        # Connect to original 2 vertices to the new centre vertex
        vc.connect(v1)
        vc.connect(v2)
        return vc

    def vpool(self, origin, supremum):
        vot = tuple(origin)
        vst = tuple(supremum)
        # Initiate vertices in case they don't exist
        vo = self.V[vot]
        vs = self.V[vst]

        # Remove origin - supremum disconnect

        # Find the lower/upper bounds of the refinement hyperrectangle
        bl = list(vot)
        bu = list(vst)
        for i, (voi, vsi) in enumerate(zip(vot, vst)):
            if bl[i] > vsi:
                bl[i] = vsi
            if bu[i] < voi:
                bu[i] = voi

        #      NOTE: This is mostly done with sets/lists because we aren't sure
        #            how well the numpy arrays will scale to thousands of
        #             dimensions.
        vn_pool = set()
        vn_pool.update(vo.nn)
        vn_pool.update(vs.nn)
        cvn_pool = copy.copy(vn_pool)
        for vn in cvn_pool:
            for i, xi in enumerate(vn.x):
                if bl[i] <= xi <= bu[i]:
                    pass
                else:
                    try:
                        vn_pool.remove(vn)
                    except KeyError:
                        pass  # NOTE: Not all neigbouds are in initial pool
        return vn_pool

    def vf_to_vv(self, vertices, simplices):
        """
        Convert a vertex-face mesh to a vertex-vertex mesh used by this class

        :param vertices: list of vertices
        :param simplices: list of simplices
        :return:
        """
        if self.dim > 1:
            for s in simplices:
                edges = itertools.combinations(s, self.dim)
                for e in edges:
                    self.V[tuple(vertices[e[0]])].connect(
                        self.V[tuple(vertices[e[1]])])
        else:
            for e in simplices:
                self.V[tuple(vertices[e[0]])].connect(
                    self.V[tuple(vertices[e[1]])])
        return

    def connect_vertex_non_symm(self, v_x, near=None):
        """
        Adds a vertex at coords v_x to the complex that is not symmetric to the
        initial triangulation and sub-triangulation.

        If near is specified (for example; a star domain or collections of
        cells known to contain v) then only those simplices containd in near
        will be searched, this greatly speeds up the process.

        If near is not specified this method will search the entire simplicial
        complex structure.

        :param v_x: tuple, coordinates of non-symmetric vertex
        :param near: set or list of vertices, these are points near v to check for
        :return:
        """
        if near is None:
            star = self.V
        else:
            star = near
        # Create the vertex origin
        if tuple(v_x) in self.V.cache:
            if self.V[v_x] in self.V_non_symm:
                pass
            else:
                return

        self.V[v_x]
        found_nn = False
        S_rows = []
        for v in star:
            S_rows.append(v.x)

        S_rows = numpy.array(S_rows)
        A = numpy.array(S_rows) - numpy.array(v_x)
        #Iterate through all the possible simplices of S_rows
        for s_i in itertools.combinations(range(S_rows.shape[0]),
                                          r=self.dim + 1):
            # Check if connected, else s_i is not a simplex
            valid_simplex = True
            for i in itertools.combinations(s_i, r=2):
                # Every combination of vertices must be connected, we check of
                # the current iteration of all combinations of s_i are connected
                # we break the loop if it is not.
                if ((not self.V[tuple(S_rows[i[1]])] in
                        self.V[tuple(S_rows[i[0]])].nn)
                    and (not (self.V[tuple(S_rows[i[0]])] in
                        self.V[tuple(S_rows[i[1]])].nn))):
                    valid_simplex = False
                    break

            S = S_rows[tuple([s_i])]
            if valid_simplex:
                if self.deg_simplex(S, proj=None):
                    valid_simplex = False

            # If s_i is a valid simplex we can test if v_x is inside si
            if valid_simplex:
                # Find the A_j0 value from the precalculated values
                A_j0 = A[[s_i]]
                if self.in_simplex(S, v_x, A_j0):
                    found_nn = True
                    break  # breaks the main for loop, s_i is the target simplex

        # Connect the simplex to point
        if found_nn:
            for i in s_i:
                self.V[v_x].connect(self.V[tuple(S_rows[i])])
        # Attached the simplex to storage for all non-symmetric vertices
        self.V_non_symm.append(self.V[v_x])
        return found_nn  # this bool value indicates a successful connection if True

    def in_simplex(self, S, v_x, A_j0=None):
        """
        Check if a vector v_x is in simplex S
        :param S: array containing simplex entries of vertices as rows
        :param v_x: a candidate vertex
        :param A_j0: array, optional, allows for A_j0 to be pre-calculated
        :return: boolean, if True v_x is in S, if False v_x is not in S

        Notes:
        https://stackoverflow.com/questions/21819132/how-do-i-check-if-a-simplex-contains-the-origin
        """
        A_11 = numpy.delete(S, 0, 0) - S[0]

        sign_det_A_11 = numpy.sign(numpy.linalg.det(A_11))
        if sign_det_A_11 == 0:
            #NOTE: We keep the variable A_11, but we loop through A_jj
            #ind=
            #while sign_det_A_11 == 0:
            #    A_11 = numpy.delete(S, ind, 0) - S[ind]
            #    sign_det_A_11 = numpy.sign(numpy.linalg.det(A_11))

            sign_det_A_11 = -1  #TODO: Choose another det of j instead?
            #TODO: Unlikely to work in many cases

        if A_j0 is None:
            A_j0 = S - v_x

        for d in range(self.dim + 1):
            det_A_jj = (-1)**d * sign_det_A_11
            #TODO: Note that scipy might be faster to add as an optional
            #      dependency
            sign_det_A_j0 = numpy.sign(numpy.linalg.det(numpy.delete(A_j0, d, 0)))
            #TODO: Note if sign_det_A_j0 == then the point is coplanar to the
            #      current simplex facet, so perhaps return True and attach?
            if det_A_jj == sign_det_A_j0:
                continue
            else:
                return False

        return True

    def deg_simplex(self, S, proj=None):
        """
        Test a simplex S for degeneracy (linear dependence in R^dim)
        :param S: Numpy array of simplex with rows as vertex vectors
        :param proj: array, optional, if the projection S[1:] - S[0] is already
                     computed it can be added as an optional argument.
        :return:
        """
        # Strategy: we test all combination of faces, if any of the determinants
        # are zero then the vectors lie on the same face and is therefore
        # linearly dependent in the space of R^dim
        if proj is None:
            proj = S[1:] - S[0]

        #TODO: Is checking the projection of one vertex against faces of other
        #       vertices sufficient? Or do we need to check more vertices in
        #       dimensions higher than 2?
        #TODO: Literature seems to suggest using proj.T, but why is this needed?
        if numpy.linalg.det(proj) == 0.0: #TODO: Repalace with tolerance?
            return True  # Simplex is degenerate
        else:
            return False  # Simplex is not degenerate



    ### TODO: DELETE ALL BELOW:
    def plot_complex(self, show=True, directed=True, complex_plot=True,
                     contour_plot=True, surface_plot=True,
                     surface_field_plot=True, minimiser_points=True,
                     point_color='do', line_color='do',
                     complex_color_f='lo', complex_color_e='do', pointsize=7,
                     no_grids=False, save_fig=True, strpath=None,
                     plot_path='fig/', fig_name='complex.pdf', arrow_width=None,
                     fig_surface=None, ax_surface=None, fig_complex=None,
                     ax_complex=None
                     ):
        """
        Plots the current simplicial complex contained in the class. It requires
        at least one vector in the self.V to have been defined.


        :param show: boolean, optional, show the output plots
        :param directed: boolean, optional, adds directed arrows to edges
        :param contour_plot: boolean, optional, contour plots of the field functions
        :param surface_plot: boolean, optional, a 3 simplicial complex + sfield plot
        :param surface_field_plot: boolean, optional, 3 dimensional surface + contour plot
        :param minimiser_points: boolean, optional, adds minimiser points
        :param point_color: str or vec, optional, colour of complex points
        :param line_color: str or vec, optional, colour of complex edges
        :param complex_color_f: str or vec, optional, colour of surface complex faces
        :param complex_color_e: str or vec, optional, colour of surface complex edges
        :param pointsize: float, optional, size of vectices on plots
        :param no_grids: boolean, optional, removes gridlines and axes
        :param save_fig: boolean, optional, save the output figure to file
        :param strpath: str, optional, string path of the file name
        :param plot_path: str, optional, relative path to file outputs
        :param fig_name: str, optional, name of the complex file to save
        :param arrow_width: float, optional, fixed size for arrows
        :return: self.ax_complex, a matplotlib Axes class containing the complex and field contour
        :return: self.ax_surface, a matplotlib Axes class containing the complex surface and field surface
        TODO: self.fig_* missing
        Examples
        --------
        # Initiate a complex class
        >>> import pylab
        >>> H = Complex(2, domain=[(0, 10)], sfield=func)

        # As an example we'll use the built in triangulation to generate vertices
        >>> H.triangulate()
        >>> H.split_generation()

        # Plot the complex
        >>> H.plot_complex()

        # You can add any sort of custom drawings to the Axes classes of the
        plots
        >>> H.ax_complex.plot(0.25, 0.25, '.', color='k', markersize=10)
        >>> H.ax_surface.scatter(0.25, 0.25, 0.25, '.', color='k', s=10)

        # Show the final plot
        >>> pylab.show()

        # Clear current plot instances
        >>> H.plot_clean()

        Example 2: Subplots  #TODO: Test
        >>> import matplotlib.pyplot as plt
        >>> fig, axes = plt.subplots(ncols=2)
        >>> H = Complex(2, domain=[(0, 10)], sfield=func)
        >>> H.triangulate()
        >>> H.split_generation()

        # Plot the complex on the same subplot
        >>> H.plot_complex(fig_surface=fig, ax_surface=axes[0],
        ...                fig_complex=fig, ax_complex=axes[1])

        # Note you can also plot several complex objects on larger subplots
        #  using this method.

        """
        if not matplotlib_available:
            logging.warning("Plotting functions are unavailable. To "
                            "install matplotlib install using ex. `pip install "
                            "matplotlib` ")
            return
        if self.sfield is None:
            directed = False  #TODO: We used this to avoid v.minimiser_point
                              # errors when is no field, should check for field
                              # instead
        # Check if fix or ax arguments are passed
        if fig_complex is not None:
            self.fig_complex = fig_complex
        if ax_complex is not None:
            self.ax_complex = ax_complex
        if fig_surface is not None:
            self.fig_surface = fig_surface
        if ax_surface is not None:
            self.ax_surface = ax_surface

        # Create pyplot.figure instance if none exists yet
        try:
            self.fig_complex
        except AttributeError:
            self.fig_complex = pyplot.figure()

        # Consistency
        if self.sfield is None:
            if contour_plot:
                contour_plot = False
                logging.warning("Warning, no associated scalar field found. "
                                "Not plotting contour_plot.")

            if surface_field_plot:
                surface_field_plot = False
                logging.warning("Warning, no associated scalar field found. "
                                "Not plotting surface field.")

        # Define colours:
        coldict = {'lo': numpy.array([242, 189, 138]) / 255,  # light orange
                   'do': numpy.array([235, 129, 27]) / 255  # Dark alert orange
                   }

        def define_cols(col):
            if (col is 'lo') or (col is 'do'):
                col = coldict[col]
            elif col is None:
                col = None
            return col

        point_color = define_cols(point_color)  # None will generate
        line_color = define_cols(line_color)
        complex_color_f = define_cols(complex_color_f)
        complex_color_e = define_cols(complex_color_e)

        if self.dim == 1:
            if arrow_width is not None:
                self.arrow_width = arrow_width
                self.mutation_scale = 58.83484054145521 * self.arrow_width * 1.3
            else:  # heuristic
                dx = self.bounds[0][1] - self.bounds[0][0]
                self.arrow_width = (dx * 0.13
                                    / (numpy.sqrt(len(self.V.cache))))
                self.mutation_scale = 58.83484054145521 * self.arrow_width * 1.3

            try:
                self.ax_complex
            except:
                self.ax_complex = self.fig_complex.add_subplot(1, 1, 1)

            min_points = []
            for v in self.V.cache:
                self.ax_complex.plot(v, 0, '.',
                                     color=point_color,
                                     markersize=pointsize)
                xlines = []
                ylines = []
                for v2 in self.V[v].nn:
                    xlines.append(v2.x)
                    ylines.append(0)

                    if directed:
                        if self.V[v].f > v2.f:  # direct V2 --> V1
                            x1_vec = list(self.V[v].x)
                            x2_vec = list(v2.x)
                            x1_vec.append(0)
                            x2_vec.append(0)
                            ap = self.plot_directed_edge(self.V[v].f, v2.f,
                                                         x1_vec, x2_vec,
                                                         mut_scale=0.5 * self.mutation_scale,
                                                         proj_dim=2,
                                                         color=line_color)

                            self.ax_complex.add_patch(ap)

                if directed:
                    if minimiser_points:
                        if self.V[v].minimiser():
                            v_min = list(v)
                            v_min.append(0)
                            min_points.append(v_min)

                self.ax_complex.plot(xlines, ylines, color=line_color)

            if minimiser_points:
                self.ax_complex = self.plot_min_points(self.ax_complex,
                                                       min_points,
                                                       proj_dim=2,
                                                       point_color=point_color,
                                                       pointsize=pointsize)

            # Clean up figure
            if self.bounds is None:
                pyplot.ylim([-1e-2, 1 + 1e-2])
                pyplot.xlim([-1e-2, 1 + 1e-2])
            else:
                fac = 1e-2  # TODO: TEST THIS
                pyplot.ylim([0 - fac, 0 + fac])
                pyplot.xlim(
                    [self.bounds[0][0] - fac * (self.bounds[0][1]
                                                - self.bounds[0][0]),
                     self.bounds[0][1] + fac * (self.bounds[0][1]
                                                - self.bounds[0][0])])
            if no_grids:
                self.ax_complex.set_xticks([])
                self.ax_complex.set_yticks([])
                self.ax_complex.axis('off')

            # Surface plots
            if surface_plot or surface_field_plot:
                try:
                    self.fig_surface
                except AttributeError:
                    self.fig_surface = pyplot.figure()
                try:
                    self.ax_surface
                except:
                    self.ax_surface = self.fig_surface.add_subplot(1, 1, 1)

                # Add a plot of the field function.
                if surface_field_plot:
                    self.fig_surface, self.ax_surface = self.plot_field_surface(
                        self.fig_surface,
                        self.ax_surface,
                        self.bounds,
                        self.sfield,
                        self.sfield_args,
                        proj_dim=2,
                        color=complex_color_f)  # TODO: Custom field colour

                if surface_plot:
                    self.fig_surface, self.ax_surface = self.plot_complex_surface(
                        self.fig_surface,
                        self.ax_surface,
                        directed=directed,
                        pointsize=pointsize,
                        color_e=complex_color_e,
                        color_f=complex_color_f,
                        min_points=min_points)

                if no_grids:
                    self.ax_surface.set_xticks([])
                    self.ax_surface.set_yticks([])
                    self.ax_surface.axis('off')

        elif self.dim == 2:
            if arrow_width is not None:
                self.arrow_width = arrow_width
            else:  # heuristic
                dx1 = self.bounds[0][1] - self.bounds[0][0]
                dx2 = self.bounds[1][1] - self.bounds[1][0]

                try:
                    self.arrow_width = (min(dx1, dx2) * 0.13
                                        / (numpy.sqrt(len(self.V.cache))))
                except TypeError:  # Allow for decimal operations
                    self.arrow_width = (min(dx1, dx2) * decimal.Decimal(0.13)
                                        / decimal.Decimal(
                                (numpy.sqrt(len(self.V.cache)))))

            try:
                self.mutation_scale = 58.8348 * self.arrow_width * 1.5
            except TypeError:  # Allow for decimal operations
                self.mutation_scale = (decimal.Decimal(58.8348)
                                       * self.arrow_width
                                       * decimal.Decimal(1.5))

            try:
                self.ax_complex
            except:
                self.ax_complex = self.fig_complex.add_subplot(1, 1, 1)

            if contour_plot:
                self.plot_contour(self.bounds, self.sfield,
                                  self.sfield_args)

            if complex_plot:
                min_points = []
                for v in self.V.cache:
                    self.ax_complex.plot(v[0], v[1], '.', color=point_color,
                                         markersize=pointsize)

                    xlines = []
                    ylines = []
                    for v2 in self.V[v].nn:
                        xlines.append(v2.x[0])
                        ylines.append(v2.x[1])
                        xlines.append(v[0])
                        ylines.append(v[1])

                        if directed:
                            if self.V[v].f > v2.f:  # direct V2 --> V1
                                ap = self.plot_directed_edge(self.V[v].f, v2.f,
                                                             self.V[v].x, v2.x,
                                                             mut_scale=self.mutation_scale,
                                                             proj_dim=2,
                                                             color=line_color)

                                self.ax_complex.add_patch(ap)
                    if directed:
                        if minimiser_points:
                            if self.V[v].minimiser():
                                min_points.append(v)

                    self.ax_complex.plot(xlines, ylines, color=line_color)

                if directed:
                    if minimiser_points:
                        self.ax_complex = self.plot_min_points(self.ax_complex,
                                                               min_points,
                                                               proj_dim=2,
                                                               point_color=point_color,
                                                               pointsize=pointsize)
                else:
                    min_points = []

            # Clean up figure
            if self.bounds is None:
                pyplot.ylim([-1e-2, 1 + 1e-2])
                pyplot.xlim([-1e-2, 1 + 1e-2])
            else:
                fac = 1e-2  # TODO: TEST THIS
                pyplot.ylim(
                    [self.bounds[1][0] - fac * (self.bounds[1][1]
                                                - self.bounds[1][0]),
                     self.bounds[1][1] + fac * (self.bounds[1][1]
                                                - self.bounds[1][0])])
                pyplot.xlim(
                    [self.bounds[0][0] - fac * (self.bounds[1][1]
                                                - self.bounds[1][0]),
                     self.bounds[0][1] + fac * (self.bounds[1][1]
                                                - self.bounds[1][0])])

            if no_grids:
                self.ax_complex.set_xticks([])
                self.ax_complex.set_yticks([])
                self.ax_complex.axis('off')

            # Surface plots
            if surface_plot or surface_field_plot:
                try:
                    self.fig_surface
                except AttributeError:
                    self.fig_surface = pyplot.figure()
                try:
                    self.ax_surface
                except:
                    self.ax_surface = self.fig_surface.gca(projection='3d')

                # Add a plot of the field function.
                if surface_field_plot:
                    self.fig_surface, self.ax_surface = self.plot_field_surface(
                        self.fig_surface,
                        self.ax_surface,
                        self.bounds,
                        self.sfield,
                        self.sfield_args,
                        proj_dim=3)

                if surface_plot:
                    self.fig_surface, self.ax_surface = self.plot_complex_surface(
                        self.fig_surface,
                        self.ax_surface,
                        directed=directed,
                        pointsize=pointsize,
                        color_e=complex_color_e,
                        color_f=complex_color_f,
                        min_points=min_points)

                if no_grids:
                    self.ax_surface.set_xticks([])
                    self.ax_surface.set_yticks([])
                    self.ax_surface.axis('off')


        elif self.dim == 3:
            try:
                self.ax_complex
            except:
                self.ax_complex = Axes3D(self.fig_complex)

            min_points = []
            for v in self.V.cache:
                self.ax_complex.scatter(v[0], v[1], v[2],
                                        color=point_color, s=pointsize)
                x = []
                y = []
                z = []
                x.append(self.V[v].x[0])
                y.append(self.V[v].x[1])
                z.append(self.V[v].x[2])
                for v2 in self.V[v].nn:
                    x.append(v2.x[0])
                    y.append(v2.x[1])
                    z.append(v2.x[2])
                    x.append(self.V[v].x[0])
                    y.append(self.V[v].x[1])
                    z.append(self.V[v].x[2])
                    if directed:
                        if self.V[v].f > v2.f:  # direct V2 --> V1
                            ap = self.plot_directed_edge(self.V[v].f, v2.f,
                                                         self.V[v].x, v2.x,
                                                         proj_dim=3,
                                                         color=line_color)
                            self.ax_complex.add_artist(ap)

                self.ax_complex.plot(x, y, z,
                                     color=line_color)
                if directed:
                    if minimiser_points:
                        if self.V[v].minimiser():
                            min_points.append(v)

            if minimiser_points:
                self.ax_complex = self.plot_min_points(self.ax_complex,
                                                       min_points,
                                                       proj_dim=3,
                                                       point_color=point_color,
                                                       pointsize=pointsize)

            self.fig_surface = None  # Current default
            self.ax_surface = None  # Current default

        else:
            logging.warning("dimension higher than 3 or wrong complex format")
            self.fig_complex = None
            self.ax_complex = None
            self.fig_surface = None
            self.ax_surface = None

        # Save figure to file
        if save_fig:
            if strpath is None:
                script_dir = os.getcwd()  # os.path.dirname(__file__)
                results_dir = os.path.join(script_dir, plot_path)
                sample_file_name = fig_name

                if not os.path.isdir(results_dir):
                    os.makedirs(results_dir)

                strpath = results_dir + sample_file_name

            self.plot_save_figure(strpath)


        if show and (not self.dim > 3):
            self.fig_complex.show()
        try:
            self.fig_surface
            self.ax_surface
            if show:
                self.fig_surface.show()
        except AttributeError:
            self.fig_surface = None  # Set to None for return reference
            self.ax_surface = None

        return self.fig_complex, self.ax_complex, self.fig_surface, self.ax_surface


    def plot_save_figure(self, strpath):

        self.fig_complex.savefig(strpath, transparent=True,
                                 bbox_inches='tight', pad_inches=0)

    def plot_clean(self, del_ax=True, del_fig=True):
        try:
            if del_ax:
                del (self.ax_complex)
            if del_fig:
                del (self.fig_complex)
        except AttributeError:
            pass

    def plot_contour(self, bounds, func, func_args=()):
        """
        Plots the field functions. Mostly for developmental purposes
        :param fig:
        :param bounds:
        :param func:
        :param func_args:
        :param surface:
        :param contour:
        :return:
        """
        xg, yg, Z = self.plot_field_grids(bounds, func, func_args)
        cs = pyplot.contour(xg, yg, Z, cmap='binary_r', color='k')
        pyplot.clabel(cs)

    def plot_complex_surface(self, fig, ax, directed=True, pointsize=5,
                             color_e=None, color_f=None, min_points=[]):
        """
        fig and ax need to be supplied outside the method
        :param fig: ex. ```fig = pyplot.figure()```
        :param ax: ex.  ```ax = fig.gca(projection='3d')```
        :param bounds:
        :param func:
        :param func_args:
        :return:
        """
        if self.dim == 1:
            # Plot edges
            z = []
            for v in self.V.cache:
                ax.plot(v, self.V[v].f, '.', color=color_e,
                        markersize=pointsize)
                z.append(self.V[v].f)
                for v2 in self.V[v].nn:
                    ax.plot([v, v2.x],
                            [self.V[v].f, v2.f],
                            color=color_e)

                    if directed:
                        if self.V[v].f > v2.f:  # direct V2 --> V1
                            x1_vec = [float(self.V[v].x[0]), self.V[v].f]
                            x2_vec = [float(v2.x[0]), v2.f]

                            a = self.plot_directed_edge(self.V[v].f, v2.f,
                                                        x1_vec, x2_vec,
                                                        proj_dim=2,
                                                        color=color_e)
                            ax.add_artist(a)

            ax.set_xlabel('$x$')
            ax.set_ylabel('$f$')

            if len(min_points) > 0:
                iter_min = min_points.copy()
                for ind, v in enumerate(iter_min):
                    min_points[ind][1] = float(self.V[v[0],].f)

                ax = self.plot_min_points(ax,
                                          min_points,
                                          proj_dim=2,
                                          point_color=color_e,
                                          pointsize=pointsize
                                          )

        elif self.dim == 2:
            # Plot edges
            z = []
            for v in self.V.cache:
                z.append(self.V[v].f)
                for v2 in self.V[v].nn:
                    ax.plot([v[0], v2.x[0]],
                            [v[1], v2.x[1]],
                            [self.V[v].f, v2.f],
                            color=color_e)

                    if directed:
                        if self.V[v].f > v2.f:  # direct V2 --> V1
                            x1_vec = list(self.V[v].x)
                            x2_vec = list(v2.x)
                            x1_vec.append(self.V[v].f)
                            x2_vec.append(v2.f)
                            a = self.plot_directed_edge(self.V[v].f, v2.f,
                                                        x1_vec, x2_vec,
                                                        proj_dim=3,
                                                        color=color_e)

                            ax.add_artist(a)

            # TODO: For some reason adding the scatterplots for minimiser spheres
            #      makes the directed edges disappear behind the field surface
            if len(min_points) > 0:
                iter_min = min_points.copy()
                for ind, v in enumerate(iter_min):
                    min_points[ind] = list(min_points[ind])
                    min_points[ind].append(self.V[v].f)

                ax = self.plot_min_points(ax,
                                          min_points,
                                          proj_dim=3,
                                          point_color=color_e,
                                          pointsize=pointsize
                                          )

            # Triangulation to plot faces
            # Compute a triangulation #NOTE: can eat memory
            self.vertex_face_mesh()

            ax.plot_trisurf(numpy.array(self.vertices_fm)[:, 0],
                            numpy.array(self.vertices_fm)[:, 1],
                            z,
                            triangles=numpy.array(self.simplices_fm_i),
                            # TODO: Select colour scheme
                            color=color_f,
                            alpha=0.4,
                            linewidth=0.2,
                            antialiased=True)

            ax.set_xlabel('$x_1$')
            ax.set_ylabel('$x_2$')
            ax.set_zlabel('$f$')

        return fig, ax

    def plot_field_surface(self, fig, ax, bounds, func, func_args=(),
                           proj_dim=2, color=None):
        """
        fig and ax need to be supplied outside the method
        :param fig: ex. ```fig = pyplot.figure()```
        :param ax: ex.  ```ax = fig.gca(projection='3d')```
        :param bounds:
        :param func:
        :param func_args:
        :return:
        """
        if proj_dim == 2:
            from matplotlib import cm
            xr = numpy.linspace(self.bounds[0][0], self.bounds[0][1], num=1000)
            fr = numpy.zeros_like(xr)
            for i in range(xr.shape[0]):
                fr[i] = func(xr[i], *func_args)

            ax.plot(xr, fr, alpha=0.6, color=color)

            ax.set_xlabel('$x$')
            ax.set_ylabel('$f$')

        if proj_dim == 3:
            from matplotlib import cm
            xg, yg, Z = self.plot_field_grids(bounds, func, func_args)
            ax.plot_surface(xg, yg, Z, rstride=1, cstride=1,
                            # cmap=cm.coolwarm,
                            # cmap=cm.magma,
                        #    cmap=cm.plasma,  #TODO: Restore
                            # cmap=cm.inferno,
                            # cmap=cm.pink,
                            # cmap=cm.viridis,
                            #facecolors="do it differently, ok?",
                            color = [0.94901961, 0.74117647, 0.54117647],
                            linewidth=0,
                            antialiased=True, alpha=0.8, shade=True)

            ax.set_xlabel('$x_1$')
            ax.set_ylabel('$x_2$')
            ax.set_zlabel('$f$')
        return fig, ax

    def plot_field_grids(self, bounds, func, func_args):
        try:
            return self.plot_xg, self.plot_yg, self.plot_Z
        except AttributeError:
            X = numpy.linspace(bounds[0][0], bounds[0][1])
            Y = numpy.linspace(bounds[1][0], bounds[1][1])
            xg, yg = numpy.meshgrid(X, Y)
            Z = numpy.zeros((xg.shape[0],
                             yg.shape[0]))

            for i in range(xg.shape[0]):
                for j in range(yg.shape[0]):
                    Z[i, j] = func(numpy.array([xg[i, j], yg[i, j]]),
                                   *func_args)

            self.plot_xg, self.plot_yg, self.plot_Z = xg, yg, Z
            return self.plot_xg, self.plot_yg, self.plot_Z

    def plot_directed_edge(self, f_v1, f_v2, x_v1, x_v2, mut_scale=20,
                           proj_dim=2,
                           color=None):
        """
        Draw a directed edge embeded in 2 or 3 dimensional space between two
        vertices v1 and v2.

        :param f_v1: field value at f(v1)
        :param f_v2: field value at f(v2)
        :param x_v1: coordinate vector 1
        :param x_v2: coordinate vector 2
        :param proj_dim: int, must be either 2 or 3
        :param color: edge color
        :return: a, artist arrow object (add with ex. Axes.add_artist(a)
        """
        if proj_dim == 2:
            if f_v1 > f_v2:  # direct V2 --> V1
                dV = numpy.array(x_v1) - numpy.array(x_v2)
                ap = matplotlib.patches.FancyArrowPatch(
                    numpy.array(x_v2) + 0.5 * dV,  # tail
                    numpy.array(x_v2) + 0.6 * dV,  # head
                    mutation_scale=mut_scale,
                    arrowstyle='-|>',
                    fc=color, ec=color,
                    color=color,
                )

        if proj_dim == 3:
            if f_v1 > f_v2:  # direct V2 --> V1
                dV = numpy.array(x_v1) - numpy.array(x_v2)
                # TODO: Might not be correct (unvalidated)
                ap = Arrow3D([x_v2[0], x_v2[0] + 0.5 * dV[0]],
                             [x_v2[1], x_v2[1] + 0.5 * dV[1]],
                             [x_v2[2], x_v2[2] + 0.5 * dV[2]],
                             mutation_scale=20,
                             lw=1, arrowstyle="-|>",
                             color=color)

        return ap

    def plot_min_points(self, axes, min_points, proj_dim=2, point_color=None,
                        pointsize=5):
        """
        Add a given list of highlighted minimiser points to axes

        :param ax: An initiated matplotlib Axes class
        :param min_points: list of minimsier points
        :param proj_dim: projection dimension, must be either 2 or 3
        :param point_color: optional
        :param point_size: optional
        :return:
        """
        if proj_dim == 2:
            for v in min_points:
                if point_color is 'r':
                    min_col = 'k'
                else:
                    min_col = 'r'

                axes.plot(v[0], v[1], '.', color=point_color,
                          markersize=2.5 * pointsize)

                axes.plot(v[0], v[1], '.', color='k',
                          markersize=1.5 * pointsize)

                axes.plot(v[0], v[1], '.', color=min_col,
                          markersize=1.4 * pointsize)

        if proj_dim == 3:
            for v in min_points:
                if point_color is 'r':
                    min_col = 'k'
                else:
                    min_col = 'r'

                axes.scatter(v[0], v[1], v[2], color=point_color,
                             s=2.5 * pointsize)

                axes.scatter(v[0], v[1], v[2], color='k',
                             s=1.5 * pointsize)

                axes.scatter(v[0], v[1], v[2], color=min_col,
                             s=1.4 * pointsize)

        return axes

    def vertex_face_mesh(self, field_conversions=True):
        """
        Convert the current simplicial complex from the default
        vertex-vertex mesh (low memory) to a

        NM

        :param field_conversions, boolean, optional
                If True then any associated field properties will be added to
                ordered lists ex, self.sfield_vf and self.vfield_vf

        :return: self.vertices_fm, A list of vertex vectors (corresponding to
                                   the ordered dict in self.V.cache)
                 self.simplices_fm, A list of (dim + 1)-lists containing vertex
                                    objects in a simplex.

                 self.simplices_fm_i, Same as self.simplices_fm except contains
                                      the indices corresponding to the list in
                                      self.vertices_fm

                 self.sfield_fm, Scalar field values corresponding to the
                                 vertices in self.vertices_fm
        """
        #TODO: UNTESindexTED FOR DIMENSIONS HIGHER THAN 2
        self.vertices_fm = []  # Vertices (A list of ob
        self.simplices_fm = []  # Faces
        self.simplices_fm_i = []

        # TODO: Add in field
        for v in self.V.cache:  # Note that cache is an OrderedDict
            self.vertices_fm.append(v)

            #TODO: Should this not be outside the initial loop?
            simplex = (self.dim + 1) * [None]  # Predetermined simplex sizes
            simplex[0] = self.V[v]

            # indexed simplices
            simplex_i = (self.dim + 1) * [None]
            simplex_i[0] = self.V[v].index

            # TODO: We need to recursively call a function that checks each nn
            #  and checks that the nn is in all parent nns (otherwise there is
            #  a deviding line in the simplex)
            # NOTE: The recursion depth is (self.dim + 1)

            # Start looping through the vertices in the star domain
            for v2 in self.V[v].nn:
                # For every v2 we loop through its neighbours in v2.nn, for
                # every v2.nn that is also in self.V[v].nn we want to build
                # simplices connecting to the current v. Note that we want to
                # try and build simplices from v, v2 and any connected
                # neighbours.
                # Since all simplices are ordered vertices, we can check that
                # all simplices are unique.
                # The reason we avoid recursion in higher dimensions is because
                # we only want to find the connecting chains v1--v2--v3--v1
                # to add to simplex stacks. Once a simplex is full we try to
                # find another chain v1--v2--vi--v1 to add to the simplex until
                # it is full, here vi is any neighbour of v2
                simplex_i[1] = v2.index
                build_simpl_i = simplex_i.copy()
                ind = 1

                if self.dim > 1:
                    for v3 in v2.nn:
                        # if v3 has a connection to v1, not in current simplex
                        # and not v1 itself:
                        if ((v3 in self.V[v].nn) and (v3 not in build_simpl_i)
                                and (v3 is not self.V[v])):
                            try:  # Fill simplex with v's neighbours until it is full
                                ind += 1
                                # (if v2.index not in build_simpl_i) and v2.index in v2.nn
                                build_simpl_i[ind] = v3.index
                            except IndexError:  # When the simplex is full
                                # ind = 1 #TODO: Check
                                # Append full simplex and create a new one
                                s_b_s_i = sorted(
                                    build_simpl_i)  # Sorted simplex indices
                                if s_b_s_i not in self.simplices_fm_i:
                                    self.simplices_fm_i.append(s_b_s_i)
                                    # TODO: Build simplices_fm
                                    # self.simplices_fm.append(s_b_s_i)

                                build_simpl_i = simplex_i.copy()
                                # Start the new simplex with current neighbour as second
                                #  entry
                                if ((v3 in self.V[v].nn) and (
                                        v3 not in build_simpl_i)
                                        and (v3 is not self.V[v])):
                                    build_simpl_i[2] = v3.index #TODO: Check missing index 1?
                                    ind = 2

                                if self.dim == 2:  # Special case, for dim > 2
                                    # it will not be full
                                    if s_b_s_i not in self.simplices_fm_i:
                                        self.simplices_fm_i.append(s_b_s_i)

                    # After loop check if we have a filled simplex
                    if len(build_simpl_i) == self.dim + 1:
                        if not (None in build_simpl_i):
                            s_b_s_i = sorted(
                                build_simpl_i)  # Sorted simplex indices
                            if s_b_s_i not in self.simplices_fm_i:
                                self.simplices_fm_i.append(s_b_s_i)

                    # NOTE: If we run out of v3 in v2.nn before a simplex is
                    # completed then there were not enough vertices to form a
                    # simplex with.

        # TODO: BUILD self.vertices_fm from  self.simplices_fm_i and
        # self.vertices_fm
        for s in self.simplices_fm_i:
            sl = []
            for i in s:
                sl.append(self.vertices_fm[i])

            # Append the newly built simplex
            self.simplices_fm.append(sl)

        return


### DELTE:
import collections, time, functools

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy as np

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)

class Arrow3D_old_broken_delete(FancyArrowPatch):
    """
    Arrow used in the plotting of 3D vecotrs

    ex.
    a = Arrow3D([0, 1], [0, 1], [0, 1], mutation_scale=20,
            lw=1, arrowstyle="-|>", color="k")
    ax.add_artist(a)
    """
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


# Note to avoid using external packages such as functools32 we use this code
# only using the standard library
def lru_cache(maxsize=255, timeout=None):
    """
    Thanks to ilialuk @ https://stackoverflow.com/users/2121105/ilialuk for
    this code snippet. Modifications by S. Endres
    """

    class LruCacheClass(object):
        def __init__(self, input_func, max_size, timeout):
            self._input_func = input_func
            self._max_size = max_size
            self._timeout = timeout

            # This will store the cache for this function,
            # format - {caller1 : [OrderedDict1, last_refresh_time1],
            #  caller2 : [OrderedDict2, last_refresh_time2]}.
            #   In case of an instance method - the caller is the instance,
            # in case called from a regular function - the caller is None.
            self._caches_dict = {}

        def cache_clear(self, caller=None):
            # Remove the cache for the caller, only if exists:
            if caller in self._caches_dict:
                del self._caches_dict[caller]
                self._caches_dict[caller] = [collections.OrderedDict(),
                                             time.time()]

        def __get__(self, obj, objtype):
            """ Called for instance methods """
            return_func = functools.partial(self._cache_wrapper, obj)
            return_func.cache_clear = functools.partial(self.cache_clear,
                                                        obj)
            # Return the wrapped function and wraps it to maintain the
            # docstring and the name of the original function:
            return functools.wraps(self._input_func)(return_func)

        def __call__(self, *args, **kwargs):
            """ Called for regular functions """
            return self._cache_wrapper(None, *args, **kwargs)

        # Set the cache_clear function in the __call__ operator:
        __call__.cache_clear = cache_clear

        def _cache_wrapper(self, caller, *args, **kwargs):
            # Create a unique key including the types (in order to
            # differentiate between 1 and '1'):
            kwargs_key = "".join(map(
                lambda x: str(x) + str(type(kwargs[x])) + str(kwargs[x]),
                sorted(kwargs)))
            key = "".join(
                map(lambda x: str(type(x)) + str(x), args)) + kwargs_key

            # Check if caller exists, if not create one:
            if caller not in self._caches_dict:
                self._caches_dict[caller] = [collections.OrderedDict(),
                                             time.time()]
            else:
                # Validate in case the refresh time has passed:
                if self._timeout is not None:
                    if (time.time() - self._caches_dict[caller][1]
                            > self._timeout):
                        self.cache_clear(caller)

            # Check if the key exists, if so - return it:
            cur_caller_cache_dict = self._caches_dict[caller][0]
            if key in cur_caller_cache_dict:
                return cur_caller_cache_dict[key]

            # Validate we didn't exceed the max_size:
            if len(cur_caller_cache_dict) >= self._max_size:
                # Delete the first item in the dict:
                try:
                    cur_caller_cache_dict.popitem(False)
                except KeyError:
                    pass
            # Call the function and store the data in the cache (call it
            # with the caller in case it's an instance function)
            if caller is not None:
                args = (caller,) + args
            cur_caller_cache_dict[key] = self._input_func(*args, **kwargs)

            return cur_caller_cache_dict[key]

    # Return the decorator wrapping the class (also wraps the instance to
    # maintain the docstring and the name of the original function):
    return (lambda input_func: functools.wraps(input_func)(
        LruCacheClass(input_func, maxsize, timeout)))


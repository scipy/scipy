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


# Main complex class:
class Complex:
    def __init__(self, dim, domain=None, sfield=None, sfield_args=(),
                 symmetry=None, constraints=None,
                 constraints_args=None, #TODO: Add
                 workers=None):
        """
        A base class for a simplicial complex described as a cache of vertices
        together with their connections.

        Important methods:
            Domain triangulation:
                    Complex.triangulate, Complex.split_generation
            Triangulating arbitrary points (must be traingulable,
                may exist outside domain):
                    Complex.triangulate(sample_set)  #TODO
            Converting another simplicial complex structure data type to the
                structure used in Complex (ex. OBJ wavefront)
                    Complex.convert(datatype, data)  #TODO
            Convert the structure in the Complex to other data type:
                    #TODO

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
            self.bounds = domain  # TODO: Assert that len(domain) is dim
        self.symmetry = symmetry  # TODO: Define the functions to be used
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

        self.H = []  # Storage structure of vertex groups #TODO: Dated?

        # Cache of all vertices
        # Initiate a vertex cache and an associated field cache, note that
        # the field case is always initiated inside the vertex cache if an
        # associated field scalar field is defined:
        if sfield is not None:
            self.V = VertexCacheField(field=sfield, field_args=sfield_args,
                                      g_cons=self.g_cons,
                                      g_cons_args=self.g_args,
                                      workers=workers)
        self.V_non_symm = []  # List of non-symmetric vertices


    def __call__(self):
        return self.H

    # %% Triangulation methods
    def cyclic_product(self, bounds, origin, supremum, centroid=True,
                       printout=False):
        #TODO: Connect every vertex with the supremum before yielding
        #TODO: Change self.origin and self. supremum to local origin
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
                #TODO: The above routine works for symm inits, so maybe test if
                # vl exists and then work from there?
                # Copy lists for iteration
                # cC0x = [x[:] for x in C0x[:i + 1]]
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

                    # Update the containers
                    C0x[i + 1].append(vu)
                    C1x[i + 1].append(a_vu)

                    # Update old containers
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
                        #    b_v.connect(ab_v)  # s-s

                            # Add new list of cross pairs
                            ab_C.append((vp[0], ab_v))

            # Printing
            if printout:
                print("=" * 19)
                print("Current symmetry group:")
                print("=" * 19)
                for v in self.V.cache:
                    self.V[v].print_out()
                print("=" * 19)

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
            # TODO: If not initial triangulation, we'll need to use a different
            # container
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

        #TODO: Add check that len(symmetry) is equal to len(self.bounds)
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
            self.cp = self.cyclic_product(cbounds, origin, supremum, centroid,
                                          printout)
            for i in self.cp:
                i

            try:
                self.triangulated_vectors.append((tuple(self.origin),
                                                  tuple(self.supremum)))
            except (AttributeError, KeyError):
                self.triangulated_vectors = [(tuple(self.origin),
                                              tuple(self.supremum))]

        else:
            #Check if generator already exists
            try:
                self.cp
            except (AttributeError, KeyError):
                self.cp = self.cyclic_product(cbounds, origin, supremum,
                                              centroid, printout)

            try:
                while len(self.V.cache) < n:
                    next(self.cp)
            except StopIteration:
                #TODO: We should maybe append and add the possibility of
                #      of starting new triangulated domains on different
                #      complexes
                try:
                    self.triangulated_vectors.append((tuple(self.origin),
                                                      tuple(self.supremum)))
                except (AttributeError, KeyError):
                    self.triangulated_vectors = [(tuple(self.origin),
                                                  tuple(self.supremum))]
        return

    def refine(self, n=1):
        #TODO: Replace with while loop checking cache size instead?
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
        #TODO: We can reduce the number of exception handling and exceptions to
        #  StopIteration if we try to generate the initial objects first and
        #  then refine only have the first few while loop runs
        while len(self.V.cache) < nt:  # while loop 1
            try:  # try 1
                # Try to access triangulated_vectors, this should only be
                # defined if an initial triangulation has already been performed
                self.triangulated_vectors
                # Try a usual iteration of the current generator, if it
                # does not exist or is exhausted then produce a new generator
                try:  # try 2
                    #next(self.rls)
                    #print(f'next(self.rls) = {next(self.rls)}')
                    next(self.rls)
                except (AttributeError, StopIteration, KeyError):
                    vp = self.triangulated_vectors[0]
                    #print(f'vp = {vp}')
                    self.rls = self.refine_local_space(*vp, bounds=self.bounds)
                    # next(self.rls)
                    #print(f'next(self.rls) = {next(self.rls)}')
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
                #self.rls =self.refine_local_space_c(*vp, vpool=vn_pool_sets[i])
                self.rls = self.refine_local_space(*vp, bounds=self.bounds)
                for i in self.rls:
                    i
        except AttributeError as ae:
            if str(ae) == "'Complex' object has no attribute " \
                          "'triangulated_vectors'":
                self.triangulate(symmetry=self.symmetry, centroid=centroids)
            else:
                raise

        return

    def refine_local_space(self, origin, supremum, bounds, centroid=1,
                              printout=0):
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

                    #TODO: Should be done in the cyclical loop somehow?
                    b_vl_c = self.split_edge(b_vu.x, b_vl.x)
                    bc_vc.connect(b_vl_c)
                    ###

                    yield b_vu
                    ba_vu = self.V[tuple(ba_vu)]
                    bc_vc.connect(ba_vu)  # Connect aN cross pairs
                    d_bc_vc.connect(ba_vu)  # Connect all to centroid

                    # Split the a + b edge of the initial triangulation:
                    os_v = self.split_edge(vectors[1].x, ba_vu.x)  # o-s
                    ss_v = self.split_edge(b_vl.x, ba_vu.x)  # s-s

                    #TODO: Should be done in the cyclical loop somehow?
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
                    #ab_C.append((bc_vc, b_vl, b_vu, ba_vl, ba_vu))
                    #ab_C.append((d_bc_vc, d_b_vl, d_b_vu, d_ba_vl, d_ba_vu))
                    ab_C.append((d_bc_vc, vectors[1], b_vl, a_vu, ba_vu))
                    ab_C.append((d_bc_vc, vl, b_vl, a_vu, ba_vu))  # = prev
                    #ab_C.append((d_bc_vc, vu, b_vu, a_vl, ba_vl))

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

                    #ab_C.append((c_vc, vl, vu, a_vl, a_vu))  # CHECK

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
                # TODO: The above routine works for symm inits, so maybe test if
                # vl exists and then work from there?
                for vectors in ab_Cc:
                    ba_vl = list(vectors[3].x)
                    ba_vu = list(vectors[4].x)
                    ba_vl[i + 1] = vut[i + 1]
                    ba_vu[i + 1] = vut[i + 1]
                    ba_vu = self.V[tuple(ba_vu)]
                    yield ba_vu
                    # Split the a + b edge of the initial triangulation:
                    d_bc_vc = self.split_edge(vectors[1].x, ba_vu.x)  # o-s
                    yield ba_vu
                    d_bc_vc.connect(vectors[1])  # Connect all to centroid
                    d_bc_vc.connect(vectors[2])  # Connect all to centroid
                    d_bc_vc.connect(vectors[3])  # Connect all to centroid
                    d_bc_vc.connect(vectors[4])  # Connect all to centroid
                    yield d_bc_vc.x
                    ba_vl = self.V[tuple(ba_vl)]
                    yield ba_vl

                    # Split remaining to centre, call this centre group "d = 0.5*a"
                    d_ba_vl = self.split_edge(vectors[3].x, ba_vl.x)
                    d_ba_vu = self.split_edge(vectors[4].x, ba_vu.x)
                    d_ba_vc = self.split_edge(d_ba_vl.x, d_ba_vu.x)
                    yield d_ba_vl
                    yield d_ba_vu
                    yield d_ba_vc
                    c_vc, vl, vu, a_vl, a_vu = vectors
                    comb = [vl, vu, a_vl, a_vu, ba_vl, ba_vu]
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
                # yield at end of loop
                c_vu = self.split_edge(vu.x, a_vu.x)
                c_vu.connect(vco)
                # Connect remaining cN group vertices
                c_vc.connect(c_vu)  # Connect cN group vertices
                yield (c_vu.x)
                Cox[i + 1].append(vu)
                Ccx[i + 1].append(c_vu)
                Cux[i + 1].append(a_vu)
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
            pass # Turn this into a logging warning?

        # Add newly triangulated vectors:
        for vs in sup_set:
            self.triangulated_vectors.append((tuple(vco.x), tuple(vs.x)))

        # Extra yield to ensure that the triangulation is completed
        if centroid:
            vcn_set = set()
            vcn_list = []
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
            print('-')
            print(f'v1.x = {v1.x}')
            vnnu = v1nn.intersection(vnn)
            print(f'vnnu = {vnnu}')

            d_v0v1 = self.split_edge(v.x, v1.x)
            for o_d_v0v1 in d_v0v1_set:
                d_v0v1.connect(o_d_v0v1)
            d_v0v1_set.add(d_v0v1)
            for v2 in vnnu:
                print(f'v2 = {v2}')
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
            #TODO: Only except float
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

        # TODO: These for loops can easily be replaced by numpy operations,
        #      tests should be run to determine which method is faster.
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

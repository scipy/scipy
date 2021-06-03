'''
Static type checking stub file for scipy/spatial/qhull.pyx
'''

from typing import List, Tuple, Any, Dict

import numpy as np
from numpy.typing import ArrayLike
from typing_extensions import final

@final
class _Qhull:
    # Read-only cython attribute that behaves, more or less, like a property
    @property
    def ndim(self) -> int: ...
    mode_option: bytes
    options: bytes
    furthest_site: bool

    def __init__(
        self,
        mode_option: bytes,
        points: np.ndarray,
        options: None | bytes = ...,
        required_options: None | bytes = ...,
        furthest_site: bool = ...,
        incremental: bool = ...,
        interior_point: None | np.ndarray = ...,
    ) -> None: ...
    def check_active(self) -> None: ...
    def close(self) -> None: ...
    def get_points(self) -> np.ndarray: ...
    def add_points(
        self,
        points: ArrayLike,
        interior_point: ArrayLike = ...
    ) -> None: ...
    def get_paraboloid_shift_scale(self) -> Tuple[float, float]: ...
    def volume_area(self) -> Tuple[float, float]: ...
    def triangulate(self) -> None: ...
    def get_simplex_facet_array(self): ...
    def get_hull_points(self) -> np.ndarray: ...
    def get_hull_facets(self) -> Tuple[List[List[int]], np.ndarray]: ...
    def get_voronoi_diagram(self) -> Tuple[
        np.ndarray,
        np.ndarray,
        List[List[int]],
        List[List[int]],
        np.ndarray,
    ]: ...
    def get_extremes_2d(self) -> np.ndarray: ...

def _get_barycentric_transforms(
    points: np.ndarray,
    simplices: np.ndarray,
    eps: float
) -> np.ndarray: ...

class _QhullUser:
    ndim: int
    npoints: int
    min_bound: np.ndarray
    max_bound: np.ndarray

    def __init__(self, qhull: _Qhull, incremental: bool = ...) -> None: ...
    def close(self) -> None: ...
    def _update(self, qhull: _Qhull) -> None: ...
    def _add_points(
        self,
        points: ArrayLike,
        restart: bool = ...,
        interior_point: ArrayLike = ...
    ) -> None: ...

class Delaunay(_QhullUser):
    furthest_site: bool
    paraboloid_scale: float
    paraboloid_shift: float
    simplices: np.ndarray
    neighbors: np.ndarray
    equations: np.ndarray
    coplanar: np.ndarray
    good: np.ndarray
    nsimplex: int
    vertices: np.ndarray

    def __init__(
        self,
        points: ArrayLike,
        furthest_site: bool = ...,
        incremental: bool = ...,
        qhull_options: None | str = ...
    ) -> None: ...
    def _update(self, qhull: _Qhull) -> None: ...
    def add_points(
        self,
        points: ArrayLike,
        restart: bool = ...
    ) -> None: ...
    @property
    def points(self) -> np.ndarray: ...
    @property
    def transform(self) -> np.ndarray: ...
    @property
    def vertex_to_simplex(self) -> np.ndarray: ...
    @property
    def vertex_neighbor_vertices(self) -> Tuple[np.ndarray, np.ndarray]: ...
    @property
    def convex_hull(self) -> np.ndarray: ...
    def find_simplex(
        self,
        xi: ArrayLike,
        bruteforce: bool = ...,
        tol: float = ...
    ) -> np.ndarray: ...
    def plane_distance(self, xi: ArrayLike) -> np.ndarray: ...
    def lift_points(self, x: ArrayLike) -> np.ndarray: ...

def tsearch(tri: Delaunay, xi: ArrayLike) -> np.ndarray: ...
def _copy_docstr(dst: object, src: object) -> None: ...

class ConvexHull(_QhullUser):
    simplices: np.ndarray
    neighbors: np.ndarray
    equations: np.ndarray
    coplanar: np.ndarray
    good: None | np.ndarray
    nsimplex: int

    def __init__(
        self,
        points: ArrayLike,
        incremental: bool = ...,
        qhull_options: None | str = ...
    ) -> None: ...
    def _update(self, qhull: _Qhull) -> None: ...
    def add_points(self, points: ArrayLike,
                   restart: bool = ...) -> None: ...
    @property
    def points(self) -> np.ndarray: ...
    @property
    def vertices(self) -> np.ndarray: ...

class Voronoi(_QhullUser):
    vertices: np.ndarray
    ridge_points: np.ndarray
    ridge_vertices: List[List[int]]
    regions: List[List[int]]
    point_region: np.ndarray
    furthest_site: bool

    def __init__(
        self,
        points: ArrayLike,
        furthest_site: bool = ...,
        incremental: bool = ...,
        qhull_options: None | str = ...
    ) -> None: ...
    def _update(self, qhull: _Qhull) -> None: ...
    def add_points(
        self,
        points: ArrayLike,
        restart: bool = ...
    ) -> None: ...
    @property
    def points(self) -> np.ndarray: ...
    @property
    def ridge_dict(self) -> Dict[Tuple[int, int], List[int]]: ...

class HalfspaceIntersection(_QhullUser):
    interior_point: np.ndarray
    dual_facets: List[List[int]]
    dual_equations: np.ndarray
    dual_points: np.ndarray
    dual_volume: float
    dual_area: float
    intersections: np.ndarray
    ndim: int
    nineq: int

    def __init__(
        self,
        halfspaces: ArrayLike,
        interior_point: ArrayLike,
        incremental: bool = ...,
        qhull_options: None | str = ...
    ) -> None: ...
    def _update(self, qhull: _Qhull) -> None: ...
    def add_halfspaces(
        self,
        halfspaces: ArrayLike,
        restart: bool = ...
    ) -> None: ...
    @property
    def halfspaces(self) -> np.ndarray: ...
    @property
    def dual_vertices(self) -> np.ndarray: ...

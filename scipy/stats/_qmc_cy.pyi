from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy.typing as npt


def _cy_wrapper_centered_discrepancy(
        sample: npt.ArrayLike, 
        iterative: bool, 
        workers: int,
) -> float: ...


def _cy_wrapper_wrap_around_discrepancy(
        sample: npt.ArrayLike,
        iterative: bool, 
        workers: int,
) -> float: ...


def _cy_wrapper_mixture_discrepancy(
        sample: npt.ArrayLike,
        iterative: bool, 
        workers: int,
) -> float: ...


def _cy_wrapper_l2_star_discrepancy(
        sample: npt.ArrayLike,
        iterative: bool,
        workers: int,
) -> float: ...


def _cy_wrapper_update_discrepancy(
        x_new_view: npt.ArrayLike,
        sample_view: npt.ArrayLike,
        initial_disc: float,
) -> float: ...

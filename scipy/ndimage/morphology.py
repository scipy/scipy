from . import _morphology


__all__ = ['iterate_structure', 'generate_binary_structure',  # noqa: F822
           'binary_erosion', 'binary_dilation', 'binary_opening',
           'binary_closing', 'binary_hit_or_miss', 'binary_propagation',
           'binary_fill_holes', 'grey_erosion', 'grey_dilation',
           'grey_opening', 'grey_closing', 'morphological_gradient',
           'morphological_laplace', 'white_tophat', 'black_tophat',
           'distance_transform_bf', 'distance_transform_cdt',
           'distance_transform_edt']


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            f"scipy.ndimage.morphology has no attribute {name}.")

    return getattr(_morphology, name)

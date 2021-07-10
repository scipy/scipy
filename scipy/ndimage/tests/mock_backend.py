import numpy as np

class _MockFunction:
    def __init__(self, return_value=None):
        self.number_calls = 0
        self.return_value = return_value
        self.last_args = ([], {})

    def __call__(self, *args, **kwargs):
        self.number_calls += 1
        self.last_args = (args, kwargs)
        return self.return_value


affine_transform = _MockFunction(np.array([[0, 0], [1, 1]]))
binary_closing = _MockFunction(np.array([[0, 0], [1, 1]]))
binary_dilation = _MockFunction(np.array([[0, 0], [1, 1]]))
binary_erosion = _MockFunction(np.array([[0, 0], [1, 1]]))
binary_fill_holes = _MockFunction(np.array([[0, 0], [1, 1]]))
binary_hit_or_miss = _MockFunction(np.array([[0, 0], [1, 1]]))
binary_propagation = _MockFunction(np.array([[0, 0], [1, 1]]))
binary_opening = _MockFunction(np.array([[0, 0], [1, 1]]))
black_tophat = _MockFunction(np.array([[0, 0], [1, 1]]))
center_of_mass = _MockFunction(np.array([[0, 0], [1, 1]]))
convolve = _MockFunction(np.array([[0, 0], [1, 1]]))
convolve1d = _MockFunction(np.array([[0, 0], [1, 1]]))
correlate = _MockFunction(np.array([[0, 0], [1, 1]]))
correlate1d = _MockFunction(np.array([[0, 0], [1, 1]]))
distance_transform_bf = _MockFunction(np.array([[0, 0], [1, 1]]))
distance_transform_cdt = _MockFunction(np.array([[0, 0], [1, 1]]))
distance_transform_edt = _MockFunction(np.array([[0, 0], [1, 1]]))
extrema = _MockFunction(np.array([[0, 0], [1, 1]]))
find_objects = _MockFunction(np.array([[0, 0], [1, 1]]))
fourier_ellipsoid = _MockFunction(np.array([[0, 0], [1, 1]]))
fourier_gaussian = _MockFunction(np.array([[0, 0], [1, 1]]))
fourier_shift = _MockFunction(np.array([[0, 0], [1, 1]]))
fourier_uniform = _MockFunction(np.array([[0, 0], [1, 1]]))
gaussian_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
gaussian_filter1d = _MockFunction(np.array([[0, 0], [1, 1]]))
gaussian_gradient_magnitude = _MockFunction(np.array([[0, 0], [1, 1]]))
gaussian_laplace = _MockFunction(np.array([[0, 0], [1, 1]]))
generate_binary_structure = _MockFunction(np.array([[0, 0], [1, 1]]))
generic_gradient_magnitude = _MockFunction(np.array([[0, 0], [1, 1]]))
generic_laplace = _MockFunction(np.array([[0, 0], [1, 1]]))
generic_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
generic_filter1d = _MockFunction(np.array([[0, 0], [1, 1]]))
geometric_transform = _MockFunction(np.array([[0, 0], [1, 1]]))
grey_closing = _MockFunction(np.array([[0, 0], [1, 1]]))
grey_dilation = _MockFunction(np.array([[0, 0], [1, 1]]))
grey_erosion = _MockFunction(np.array([[0, 0], [1, 1]]))
grey_opening = _MockFunction(np.array([[0, 0], [1, 1]]))
histogram = _MockFunction(np.array([[0, 0], [1, 1]]))
iterate_structure = _MockFunction(np.array([[0, 0], [1, 1]]))
label = _MockFunction(np.array([[0, 0], [1, 1]]))
labeled_comprehension = _MockFunction(np.array([[0, 0], [1, 1]]))
laplace = _MockFunction(np.array([[0, 0], [1, 1]]))
map_coordinates = _MockFunction(np.array([[0, 0], [1, 1]]))
maximum = _MockFunction(np.array([[0, 0], [1, 1]]))
maximum_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
maximum_filter1d = _MockFunction(np.array([[0, 0], [1, 1]]))
maximum_position = _MockFunction(np.array([[0, 0], [1, 1]]))
mean = _MockFunction(np.array([[0, 0], [1, 1]]))
median = _MockFunction(np.array([[0, 0], [1, 1]]))
median_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
minimum = _MockFunction(np.array([[0, 0], [1, 1]]))
minimum_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
minimum_filter1d = _MockFunction(np.array([[0, 0], [1, 1]]))
minimum_position = _MockFunction(np.array([[0, 0], [1, 1]]))
morphological_gradient = _MockFunction(np.array([[0, 0], [1, 1]]))
morphological_laplace = _MockFunction(np.array([[0, 0], [1, 1]]))
percentile_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
prewitt = _MockFunction(np.array([[0, 0], [1, 1]]))
rank_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
rotate = _MockFunction(np.array([[0, 0], [1, 1]]))
shift = _MockFunction(np.array([[0, 0], [1, 1]]))
sobel = _MockFunction(np.array([[0, 0], [1, 1]]))
spline_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
spline_filter1d = _MockFunction(np.array([[0, 0], [1, 1]]))
standard_deviation = _MockFunction(np.array([[0, 0], [1, 1]]))
sum_labels = _MockFunction(np.array([[0, 0], [1, 1]]))
uniform_filter = _MockFunction(np.array([[0, 0], [1, 1]]))
uniform_filter1d = _MockFunction(np.array([[0, 0], [1, 1]]))
variance = _MockFunction(np.array([[0, 0], [1, 1]]))
watershed_ift = _MockFunction(np.array([[0, 0], [1, 1]]))
white_tophat = _MockFunction(np.array([[0, 0], [1, 1]]))
zoom = _MockFunction(np.array([[0, 0], [1, 1]]))


__ua_domain__ = "numpy.scipy.ndimage"


def __ua_function__(method, args, kwargs):
    fn = globals().get(method.__name__)
    return (fn(*args, **kwargs) if fn is not None
            else NotImplemented)

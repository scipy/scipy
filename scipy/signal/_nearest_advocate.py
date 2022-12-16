import numpy as np
from ._nearest_advocate_util import _nearest_advocate

def nearest_advocate(arr_ref, arr_sig,
                     td_min, td_max, sps=10, sparse_factor=1,
                     dist_max=0.0, regulate_paddings=True, dist_padding=0.0):
    """Post-hoc synchronization method for event-based time-series data.
    
    Calculates the synchronicity of two arrays of timestamps for a search space between td_min and td_max with a precision of 1/sps. The synchronicity is given by the mean of all minimal distances between each event in arr_sig and its nearest advocate in arr_ref.

    Parameters
    ----------
    arr_ref : array_like
        Sorted reference array (1-D) with timestamps assumed to be correct.
    arr_sig : array_like
        Sorted signal array (1-D) of timestamps, assumed to be shifted by an unknown constant time-delta.    
    td_min : float
        Lower bound of the search space for the time-shift.
    td_max : float 
        Upper bound of the search space for the time-shift.
    sps : int, optional
        Number of investigated time-shifts per second, should be higher than 10 times the number of median gap of `arr_ref` (default 10).
    sparse_factor : int, optional
        Factor for the sparseness of `arr_sig` for the calculation, higher is faster but may be less accurate (default 1).
    dist_max : float, optional
        Maximal accepted distances between two advocate events. It should be around 1/4 of the median gap of `arr_ref` (default).
    regulate_paddings : bool, optional
        Regulate non-overlapping events in `arr_sig` with a maximum distance of dist_padding (default True).
    dist_padding : float, optional
        Distance assigned to non-overlapping (padding) events. It should be around 1/4 of the median gap of `arr_ref` (default). Obsolete if `regulate_paddings` is False

    Returns
    -------
    time_shifts : array_like
        Two-columned 2-D array with evaluated time-shifts (between `td_min` and `td_max`) and the respective mean distances. The time-delta with the lowest mean distance is the estimation for the time-shift between the two arrays.

    Notes
    -----
    .. versionadded:: 1.9.4

    References
    ----------
    C. Schranz, S. Mayr, "Ein neuer Algorithmus zur Zeitsynchronisierung von Ereignis-basierten Zeitreihendaten als Alternative zur Kreuzkorrelation", Spinfortec (Chemnitz 2022). :doi:`10.5281/zenodo.7370958`

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(0)
    >>> from scipy.signal import nearest_advocate
    >>> N = 10_000 
    
    Create a reference array whose events differences are sampled from a normal distribution. The signal array is the reference but shifted by `np.pi` and addional gaussian noise. The event-timestamps of both arrays must be sorted.
    
    >>> arr_ref = np.sort(np.cumsum(np.random.normal(loc=1, scale=0.25, size=N))).astype(np.float32)
    >>> arr_sig = np.sort(arr_ref + np.pi + np.random.normal(loc=0, scale=0.1, size=N)).astype(np.float32)

    The function `nearest_advocate` returns a two-columned array with all investigated time-shifts and their mean distances, i.e., the measure of the synchronicity between both array (lower is better). 
    
    >>> time_shifts = nearest_advocate(arr_ref=arr_ref, arr_sig=arr_sig, td_min=-60, td_max=60, sps=20)
    >>> time_shift, min_mean_dist = time_shifts[np.argmin(time_shifts[:,1])]
    >>> print(time_shift, min_mean_dist)
    3.15, 0.079508
    
    Plot the resulting table
    
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(time_shifts[:,0], time_shifts[:,1], color="steelblue", label="Mean distance")
    >>> plt.vlines(x=time_shift, ymin=0.05, ymax=0.25, color="firebrick", label=f"Shift = {time_shift:.2f}s")
    >>> plt.xlim(0, 8)
    >>> plt.xlabel("Time shift (s)")
    >>> plt.ylabel("Mean distance (s)")
    >>> plt.legend(loc="lower right")
    >>> plt.show()
    """
    assert isinstance(arr_ref, np.ndarray) and len(arr_ref.shape)==1
    assert isinstance(arr_sig, np.ndarray) and len(arr_sig.shape)==1
    assert isinstance(td_min, (int, float))
    assert isinstance(td_max, (int, float))
    assert isinstance(sps, int)
    assert isinstance(sparse_factor, int)
    assert isinstance(dist_max, float)
    assert isinstance(regulate_paddings, bool)
    assert isinstance(dist_padding, float)
    
    # call and return the cython function
    return _nearest_advocate(arr_ref, arr_sig, 
                             td_min=td_min, td_max=td_max, sps=sps, sparse_factor=sparse_factor, 
                             dist_max=dist_max, regulate_paddings=regulate_paddings, 
                             dist_padding=dist_padding)


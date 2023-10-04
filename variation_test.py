import numpy as np
from scipy import stats
from numpy.ma import masked_array as ma

fun = stats.variation

def try_case(x, message, kwds):
    try:
        print(message)
        res = fun(x, **kwds)
        print(res)
        print(type(res), res.dtype)
    except AttributeError:
        print(type(res))
    except Exception as e:
        print(f"{type(e)}: {str(e)}")

    try:
        print("vs (master)")
        res = fun(x, **kwds,_no_deco=True)
        print(res)
        print(type(res), res.dtype)
    except AttributeError:
        print(type(res))
    except Exception as e:
        print(f"{type(e)}: {str(e)}")

    print("----")

np.random.seed(0)
random_array = np.random.rand(4, 6)
random_mask = np.random.rand(4, 6) > 0.8
random_nan_array = random_array.copy()
random_nan_array[random_mask] = np.nan
random_ma = np.ma.masked_array(random_array, mask=random_mask)
random_ma32 = np.ma.masked_array(random_array, mask=random_mask, dtype=np.float32)

cases = [
          ([], 'Zero Observations:', {}),
          ([1.], 'One Observation:', {}),
          ([1., 2.], 'Two Observations', {}),
          (random_array, '2D Array, no axis', {}),
          (random_array, '2D Array, axis=0', {'axis': 0}),
          (random_array, '2D Array, axis=1', {'axis': 1}),
          (np.matrix(random_array), 'Matrix, axis=1', {'axis': 1}, ),
          (random_array, '2D Array, axis=None', {'axis': None}),
          (random_array, '2D Array, axis=(0,1)', {'axis': (0, 1)}),
          ([np.nan, 1., 2., 3.], '1D Array with NaN, No nan_policy', {}),
          ([np.nan, 1., 2., 3.], '1D Array with NaN, raise', {'nan_policy': 'raise'}),
          ([np.nan, 1., 2., 3.], '1D Array with NaN, propagate', {'nan_policy': 'propagate'}),
          ([np.nan, 1., 2., 3.], '1D Array with NaN, omit', {'nan_policy': 'omit'}),
          ([np.nan, np.nan, np.nan], '1D Array with All NaNs, omit', {'nan_policy': 'omit'}),
          ([np.nan, np.nan, 3.], '1D Array with all but one NaN, omit', {'nan_policy': 'omit'}),
          (ma([0., 1., 2., 3.], mask=[True, False, False, False]), '1D Masked array, one masked', {}),
          (ma([1., 2., 3.], mask=[True, True, True]), '1D Masked array, all masked', {}),
          (ma([1., 2., 3.], mask=[True, True, False]), '1D Masked array, all but one masked', {}),
          (random_nan_array, '2D Array with NaNs, propagate, axis=1', {'axis': 1, 'nan_policy': 'propagate'}, ),
          (random_nan_array, '2D Array with NaNs, omit, axis=1', {'axis': 1, 'nan_policy': 'omit'}, ),
          (random_ma, '2D Masked array, axis=1', {'axis': 1}, ),
          (random_array.astype(np.float32), '2D Array, float32, axis=1', {'axis': 1}),
          (random_nan_array.astype(np.float32), '2D Array with NaNs, float32, propagate, axis=1', {'axis': 1, 'nan_policy': 'propagate'}, ),
          (random_ma32, '2D Masked Array, float32, axis=1', {'axis': 1}, ),
          ]

for case in cases:
    try_case(*case)

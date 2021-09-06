.. _uarray:

``__ua_convert__`` in ``uarray`` backend
=========================================

SciPy supports ``uarray`` in `scipy.ndimage` and `scipy.fft`.
Any ``uarray`` backend needs to define three protocols: ``__ua_domain__``,
``__ua_function__`` and ``__ua_convert__``. The first two are mandatory
and the last is optional. In this section, we provide SciPy developer
guidelines for how the optional protocol ``__ua_convert__`` is defined
and intended to be used, to ensure that as this protocol is added to other
libraries making use of ``uarray``, we maintain a consistent API.


``__ua_convert__`` is optional, but when adding ``uarray`` support in
`scipy.ndimage` backend, we need to support and accept some arguments
which can come in multiple ``types``. A couple of these common arguments
in the ndimage module functions are : ``'output'`` and ``'index'``.

* ``output``:
  This is the output array which holds the output of a function.
  It can take different types include a python ``str`` when output='f',
  or ``output='float32'`` and dtypes like ``output=np.float32``.
  Hence we need to special case these different types in ``__ua_convert__``.
  Here we create a class called ``ndimage_output`` to represent the
  output types and then check for a particular dispatch type within the
  protocol::

      if dispatch_type is ndimage_output:
          if isinstance(value, np.ndarray):
              return np.asarray(value)
          else:
              try:
                  return np.dtype(str(value))
              except TypeError:
                  return np.dtype(value)


* ``index``:
  This argument in most of the ``ndimage.measurements`` submodule functions
  accepts an integer as well as an an array. Similar to ``ndimage_output``
  we have an ``ndimage_index`` class to handle both different types which
  the ``index`` argument supports, and then check for a particular dispatch
  type within the protocol::

      if dispatch_type is ndimage_index:
          if isinstance(value, np.ndarray):
              return np.asarray(value)
          else:
              return value

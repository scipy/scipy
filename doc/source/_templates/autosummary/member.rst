:orphan:

{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

member

.. auto{{ objtype }}:: {{ fullname | replace("scipy.", "scipy::") }}

{# In the fullname (e.g. `numpy.ma.MaskedArray.methodname`), the module name
is ambiguous. Using a `::` separator (e.g. `numpy::ma.MaskedArray.methodname`)
specifies `numpy` as the module name. #}

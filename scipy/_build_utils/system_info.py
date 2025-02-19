def combine_dict(*dicts, **kw):
    """
    Combine Numpy distutils style library configuration dictionaries.

    Parameters
    ----------
    *dicts
        Dictionaries of keys. List-valued keys will be concatenated.
        Otherwise, duplicate keys with different values result to
        an error. The input arguments are not modified.
    **kw
        Keyword arguments are treated as an additional dictionary
        (the first one, i.e., prepended).

    Returns
    -------
    combined
        Dictionary with combined values.
    """
    new_dict = {}

    for d in (kw,) + dicts:
        for key, value in d.items():
            if new_dict.get(key, None) is not None:
                old_value = new_dict[key]
                if isinstance(value, list | tuple):
                    if isinstance(old_value, list | tuple):
                        new_dict[key] = list(old_value) + list(value)
                        continue
                elif value == old_value:
                    continue

                raise ValueError(f"Conflicting configuration dicts: {new_dict!r} {d!r}")
            else:
                new_dict[key] = value

    return new_dict

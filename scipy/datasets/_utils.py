import os
import shutil
from ._registry import registry

try:
    import appdirs
except ImportError:
    appdirs = None


def _clear_cache(datasets, cache_dir=None, data_registry=None):
    if data_registry is None:
        # Use SciPy Datasets registry map
        data_registry = registry
    if cache_dir is None:
        # Use default cache_dir path
        if appdirs is None:
            # appdirs is pooch dependency
            raise ImportError("Missing optional dependency 'pooch' required "
                              "for scipy.datasets module. Please use pip or "
                              "conda to install 'pooch'.")
        cache_dir = appdirs.user_cache_dir("scipy-data")

    if not os.path.exists(cache_dir):
        print(f"Cache Directory {cache_dir} doesn't exist. Nothing to clear.")
        return

    if datasets is None:
        print(f"Cleaning the cache directory {cache_dir}!")
        shutil.rmtree(cache_dir)
    else:
        for dataset_name in datasets:
            if dataset_name not in data_registry:
                raise ValueError(f"Dataset {dataset_name} doesn't exist. "
                                 "Please check if the passed dataset "
                                 "is a subset of the following: "
                                 f"{list(data_registry.keys())}")

            dataset_path = os.path.join(cache_dir, dataset_name)
            if os.path.exists(dataset_path):
                print(f"Cleaning the {dataset_name} file!")
                os.remove(dataset_path)
            else:
                print(f"Path {dataset_path} doesn't exist. Nothing to clear.")


def clear_cache(datasets=None):
    """
    Cleans the scipy datasets cache directory.

    If a list/tuple of dataset strings is provided, then it
    removes all the datasets in the list/tuple.

    By default, it removes all the cached data files.

    Parameters
    ----------
    datasets : list/tuple of str or None
        A list/tuple of datasets to be removed from cache.

    """
    _clear_cache(datasets)

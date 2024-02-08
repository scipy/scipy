import os
import shutil
from ._registry import method_files_map, _CACHE_DIR_FILE

try:
    import platformdirs
except ImportError:
    platformdirs = None  # type: ignore[assignment]


class CacheDir:
    """
    When a custom cache directory is required for scipy.datasets,
    an object of CacheDir may be used to read the currently set cache
    path or set it to a custom user defined path.

    Examples
    --------
    >>> from scipy import datasets
    >>> cache_dir = datasets.CacheDir()
    >>> cache_path = cache_dir.get_cache_path()
    >>> datasets.CacheDir().set_cache_path('./data_tmp')
    >>> default_cache_path = cache_dir.default_cache_path
    """
    def __init__(self):
        if platformdirs is None:
            # platformdirs is pooch dependency
            raise ImportError("Missing optional dependency 'pooch' required "
                              "for scipy.datasets module. Please use pip or "
                              "conda to install 'pooch'.")
        # Use the default cache folder for the operating system. Pooch uses
        # platformdirs (https://github.com/platformdirs/platformdirs)
        # to select an appropriate directory for the cache on each platform.
        self.default_cache_path = platformdirs.user_cache_dir("scipy-data")
        self._cache_path = self.default_cache_path

    def get_cache_path(self):
        if os.path.exists(_CACHE_DIR_FILE):
            with open(_CACHE_DIR_FILE) as f:
                self._cache_path = f.read().strip()
            return self._cache_path
        else:
            return self.default_cache_path

    def set_cache_path(self, new_path=None):
        if new_path:
            self._cache_path = os.path.abspath(new_path)
        else:
            # fallback to default cache path
            self._cache_path = self.default_cache_path
        print(f"Cache path set to {self._cache_path}")
        with open(_CACHE_DIR_FILE, 'w') as f:
            f.write(self._cache_path+'\n')


def _clear_cache(datasets, cache_dir=None, method_map=None):
    if method_map is None:
        # Use SciPy Datasets method map
        method_map = method_files_map
    if cache_dir is None:
        cache_dir = CacheDir().get_cache_path()

    if not os.path.exists(cache_dir):
        print(f"Cache Directory {cache_dir} doesn't exist. Nothing to clear.")
        return

    if datasets is None:
        print(f"Cleaning the cache directory {cache_dir}!")
        shutil.rmtree(cache_dir)
        if os.path.exists(_CACHE_DIR_FILE):
            # Remove temporary file containing cache path
            print(f"Cleaning the tmp file {_CACHE_DIR_FILE}!")
            os.remove(_CACHE_DIR_FILE)
    else:
        if not isinstance(datasets, (list, tuple)):
            # single dataset method passed should be converted to list
            datasets = [datasets, ]
        for dataset in datasets:
            assert callable(dataset)
            dataset_name = dataset.__name__  # Name of the dataset method
            if dataset_name not in method_map:
                raise ValueError(f"Dataset method {dataset_name} doesn't "
                                 "exist. Please check if the passed dataset "
                                 "is a subset of the following dataset "
                                 f"methods: {list(method_map.keys())}")

            data_files = method_map[dataset_name]
            data_filepaths = [os.path.join(cache_dir, file)
                              for file in data_files]
            for data_filepath in data_filepaths:
                if os.path.exists(data_filepath):
                    print("Cleaning the file "
                          f"{os.path.split(data_filepath)[1]} "
                          f"for dataset {dataset_name}")
                    os.remove(data_filepath)
                else:
                    print(f"Path {data_filepath} doesn't exist. "
                          "Nothing to clear.")


def clear_cache(datasets=None):
    """
    Cleans the scipy datasets cache directory.

    If a scipy.datasets method or a list/tuple of the same is
    provided, then clear_cache removes all the data files
    associated to the passed dataset method callable(s).

    By default, it removes all the cached data files.

    Parameters
    ----------
    datasets : callable or list/tuple of callable or None

    Examples
    --------
    >>> from scipy import datasets
    >>> ascent_array = datasets.ascent()
    >>> ascent_array.shape
    (512, 512)
    >>> datasets.clear_cache([datasets.ascent])
    Cleaning the file ascent.dat for dataset ascent
    """
    _clear_cache(datasets)

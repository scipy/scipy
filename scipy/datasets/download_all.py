#!python3
"""
Platform independent script to download all the
`scipy.dataset` module data files.
This doesn't require a full scipy build.

Run: python download_all.py <download_dir>
"""

import argparse
import pooch

if __package__ is None or __package__ == '':
    # Running as python script, use absolute import
    import _registry  # type: ignore
else:
    # Running as python module, use relative import
    from . import _registry


def download_all(path=pooch.os_cache('scipy-data')):
    """
    Utility method to download all the dataset files
    for `scipy.datasets` module.

    Parameters
    ----------
    path : str
        Directory path to download all the dataset files.
        Defaults to the system cache_dir detected by pooch.
    """
    for dataset_name, dataset_hash in _registry.registry.items():
        pooch.retrieve(url=_registry.registry_urls[dataset_name],
                       known_hash=dataset_hash,
                       fname=dataset_name, path=path)


def main():
    parser = argparse.ArgumentParser(description='Download SciPy dataset files.')
    parser.add_argument("path", type=str, default=pooch.os_cache('scipy-data'),
                        help="Directory path to download all the dataset files.")
    args = parser.parse_args()
    download_all(args.path)


if __name__ == "__main__":
    main()

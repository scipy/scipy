#!/usr/bin/env python3
"""
Download wheels from Anaconda staging area.
"""

import os
import re
import shutil
import argparse

import urllib3
from bs4 import BeautifulSoup

__version__ = '0.2'


def get_wheel_names(version, staging_url, prefix):
    """ Get wheel names from Anaconda HTML directory.

    This looks in the Anaconda multibuild-wheels-staging page and
    parses the HTML to get all the wheel names for a release version.

    Parameters
    ----------
    version : str
        The release version. For instance, "1.5.0".
    staging_url : str
        URL at which to find packages
    prefix : str
        Prefix for wheels to download.
    """
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED')
    tmpl = re.compile(rf"^.*{prefix}-{version}-.*\.whl$")
    index_url = f"{staging_url}/files"
    index_html = http.request('GET', index_url)
    soup = BeautifulSoup(index_html.data, 'html.parser')
    return soup.findAll(text=tmpl)


def download_wheels(version, wheelhouse, staging_url, prefix):
    """Download release wheels.

    The release wheels for the given package version are downloaded
    into the given directory.

    Parameters
    ----------
    version : str
        The release version. For instance, "1.5.0".
    wheelhouse : str
        Directory in which to download the wheels.
    staging_url : str
        URL at which to find packages
    prefix : str
        Prefix for wheels to download.
    """
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED')
    wheel_names = get_wheel_names(version, staging_url, prefix)

    for i, wheel_name in enumerate(wheel_names):
        wheel_url = f"{staging_url}/{version}/download/{wheel_name}"
        wheel_path = os.path.join(wheelhouse, wheel_name)
        with open(wheel_path, 'wb') as f:
            with http.request('GET', wheel_url, preload_content=False,) as r:
                print(f"{i + 1:<4}{wheel_name}")
                shutil.copyfileobj(r, f)
    print(f"\nTotal files downloaded: {len(wheel_names)}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "version",
        help="Package version to download.")
    parser.add_argument(
        "--staging-url",
        default='https://anaconda.org/multibuild-wheels-staging/scipy',
        help="URL at which to find packages")
    parser.add_argument(
        "--prefix",
        default='scipy',
        help="Prefix for wheels (e.g 'scipy'")
    parser.add_argument(
        "-w", "--wheelhouse",
        default=os.path.join(os.getcwd(), "release", "installers"),
        help="Directory in which to store downloaded wheels\n"
             "[defaults to <cwd>/release/installers]")

    args = parser.parse_args()

    wheelhouse = os.path.expanduser(args.wheelhouse)
    if not os.path.isdir(wheelhouse):
        raise RuntimeError(
            f"{wheelhouse} wheelhouse directory is not present."
            " Perhaps you need to use the '-w' flag to specify one.")

    download_wheels(args.version, wheelhouse, args.staging_url, args.prefix)

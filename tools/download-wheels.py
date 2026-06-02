#!/usr/bin/env python
"""
Download SciPy wheels from Anaconda staging area.

"""
import os
import re
import shutil
import argparse
import urllib
import urllib.request

import urllib3
from bs4 import BeautifulSoup

__version__ = '0.1'

# Edit these for other projects.
STAGING_FILE_URL = "https://pypi.anaconda.org/multibuild-wheels-staging/simple/scipy/"
STAGING_URL = 'https://anaconda.org/multibuild-wheels-staging/scipy'
PREFIX = 'scipy'

def http_manager():
    """
    Return a urllib3 http request manager, leveraging
    proxy settings when available.
    """
    proxy_dict = urllib.request.getproxies()
    if 'http' in proxy_dict:
        http = urllib3.ProxyManager(proxy_dict['http'])
    elif 'all' in proxy_dict:
        http = urllib3.ProxyManager(proxy_dict['all'])
    else:
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED')
    return http

def get_wheel_names(version):
    """ Get wheel names from Anaconda HTML directory.

    This looks in the Anaconda multibuild-wheels-staging page and
    parses the HTML to get all the wheel names for a release version.

    Parameters
    ----------
    version : str
        The release version. For instance, "1.5.0".

    """
    http = http_manager()
    tmpl = re.compile(rf"^.*{PREFIX}-{version}-.*\.whl$")
    index_url = f"{STAGING_FILE_URL}"
    index_html = http.request('GET', index_url)
    soup = BeautifulSoup(index_html.data, 'html.parser')
    return soup.findAll(string=tmpl)


def download_wheels(version, wheelhouse):
    """Download release wheels.

    The release wheels for the given SciPy version are downloaded
    into the given directory.

    Parameters
    ----------
    version : str
        The release version. For instance, "1.5.0".
    wheelhouse : str
        Directory in which to download the wheels.

    """
    http = http_manager()
    wheel_names = get_wheel_names(version)

    for i, wheel_name in enumerate(wheel_names):
        wheel_url = f"{STAGING_URL}/{version}/download/{wheel_name}"
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
        help="SciPy version to download.")
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

    download_wheels(args.version, wheelhouse)

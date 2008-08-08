from setuptools import setup, find_packages

setup(
    author       = 'Ilan Schnell',
    author_email = 'ischnell@enthought.com',
    description  = 'C compiled UFuncs from python source',

    name         = "mkufunc",
    version      = "0.1",

    zip_safe = False,
    package_data = {'': ['*.h']},
    packages = find_packages(),
    install_requires = ['scipy >= 0.6.0']
    )

import os
import importlib


def configure_scipy_openblas(blas_variant='32'):
    """Create .openblas/scipy-openblas.pc

    Requires a pre-installed scipy-openblas32 wheel from PyPI.
    """
    basedir = os.getcwd()
    openblas_dir = os.path.join(basedir, ".openblas")
    pkg_config_fname = os.path.join(openblas_dir, "scipy-openblas.pc")

    if os.path.exists(pkg_config_fname):
        return None

    module_name = f"scipy_openblas{blas_variant}"
    try:
        openblas = importlib.import_module(module_name)
    except ModuleNotFoundError:
        raise RuntimeError(f"'pip install {module_name} first")

    os.makedirs(openblas_dir, exist_ok=True)
    with open(pkg_config_fname, "wt", encoding="utf8") as fid:
        fid.write(openblas.get_pkg_config().replace("\\", "/"))

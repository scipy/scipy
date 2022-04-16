import os


def unuran_pre_build_hook(build_clib, build_info):
    from scipy._build_utils.compiler_helper import (get_c_std_flag,
                                                    try_compile, has_flag)
    c = build_clib.compiler
    c_flag = get_c_std_flag(c)
    if c_flag is not None:
        if "extra_compiler_args" not in build_info:
            build_info["extra_compiler_args"] = []
        build_info["extra_compiler_args"].append(c_flag)
    deps = {"unistd.h": ["HAVE_DECL_GETOPT", "HAVE_UNISTD_H"],
            "dlfcn.h": ["HAVE_DLFCN_H"],
            "sys/time.h": ["HAVE_GETTIMEOFDAY", "HAVE_SYS_TIME_H",
                           "TIME_WITH_SYS_TIME"],
            "memory.h": ["HAVE_MEMORY_H"],
            "strings.h": ["HAVE_STRCASECMP", "HAVE_STRINGS_H"],
            "sys/stat.h": ["HAVE_SYS_STAT_H"],
            "sys/types.h": ["HAVE_SYS_TYPES_H"]}
    for dep in deps:
        has_dep = try_compile(c, code=f"#include <{dep}>\n"
                                      "int main(int argc, char **argv){}")
        if has_dep:
            for macro in deps[dep]:
                build_info["macros"].append((macro, "1"))
    if has_flag(c, flag="-lm"):
        try:
            build_info["libraries"].append("m")
        except KeyError:
            build_info["libraries"] = ["m"]


def _get_sources(dirs):
    sources = []
    for dir_ in dirs:
        files = [
            file for file in os.listdir(dir_) if (not os.path.isdir(file))
        ]
        path = [str(dir_ / file) for file in files]
        sources += [source for source in path if (source.endswith(".c"))]
    return sources


def _get_version(unuran_dir, configure_dot_ac, target_name):
    configure_dot_ac = unuran_dir / configure_dot_ac
    with open(configure_dot_ac, "r") as f:
        s = f.read()
        start_idx = s.find(target_name)
        end_idx = s[start_idx:].find(")") + len(s[:start_idx])
        version = s[start_idx:end_idx].split(",")[1][1:-1]
    return version


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration
    from scipy._lib._unuran_utils import _unuran_dir

    if not os.path.exists(_unuran_dir(ret_path=True) / 'README.md'):
        raise RuntimeError("Missing the `unuran` submodule! Run `git "
                           "submodule update --init` to fix this.")

    config = Configuration("_unuran", parent_package, top_path)

    # UNU.RAN info
    UNURAN_DIR = _unuran_dir(ret_path=True).resolve()
    UNURAN_VERSION = _get_version(UNURAN_DIR, "unuran/configure.ac",
                                  "AM_INIT_AUTOMAKE")

    DEFINE_MACROS = [
        ("HAVE_ALARM", "1"),
        ("HAVE_DECL_ALARM", "1"),
        ("HAVE_DECL_HUGE_VAL", "1"),
        ("HAVE_DECL_INFINITY", "1"),
        ("HAVE_DECL_ISFINITE", "0"),
        ("HAVE_DECL_ISINF", "0"),
        ("HAVE_DECL_ISNAN", "1"),
        ("HAVE_DECL_LOG1P", "1"),
        ("HAVE_DECL_SIGNAL", "1"),
        ("HAVE_DECL_SNPRINTF", "1"),
        ("HAVE_DECL_VSNPRINTF", "1"),
        ("HAVE_FLOAT_H", "1"),
        ("HAVE_FLOOR", "1"),
        ("HAVE_IEEE_COMPARISONS", "1"),
        ("HAVE_INTTYPES_H", "1"),
        ("HAVE_LIBM", "1"),
        ("HAVE_LIMITS_H", "1"),
        ("HAVE_POW", "1"),
        ("HAVE_SIGNAL", "1"),
        ("HAVE_SQRT", "1"),
        ("HAVE_STDINT_H", "1"),
        ("HAVE_STDLIB_H", "1"),
        ("HAVE_STRCHR", "1"),
        ("HAVE_STRING_H", "1"),
        ("HAVE_STRTOL", "1"),
        ("HAVE_STRTOUL", "1"),
        ("LT_OBJDIR", '".libs/"'),
        ("PACKAGE", '"unuran"'),
        ("PACKAGE_BUGREPORT", '"unuran@statmath.wu.ac.at"'),
        ("PACKAGE_NAME", '"unuran"'),
        ("PACKAGE_STRING", '"unuran %s"' % UNURAN_VERSION),
        ("PACKAGE_TARNAME", '"unuran"'),
        ("PACKAGE_URL", '""'),
        ("PACKAGE_VERSION", '"%s"' % UNURAN_VERSION),
        ("STDC_HEADERS", "1"),
        ("UNUR_ENABLE_INFO", "1"),
        ("VERSION", '"%s"' % UNURAN_VERSION),
        ("HAVE_CONFIG_H", "1"),
        ("_ISOC99_SOURCE", "1"),
    ]

    UNURAN_DIRS = [
        os.path.join("unuran", "src"),
        os.path.join("unuran", "src", "distr"),
        os.path.join("unuran", "src", "distributions"),
        os.path.join("unuran", "src", "methods"),
        os.path.join("unuran", "src", "parser"),
        os.path.join("unuran", "src", "specfunct"),
        os.path.join("unuran", "src", "urng"),
        os.path.join("unuran", "src", "utils"),
        os.path.join("unuran", "src", "tests"),
    ]
    UNURAN_SOURCE_DIRS = [UNURAN_DIR / dir_ for dir_ in UNURAN_DIRS]

    sources = _get_sources(UNURAN_SOURCE_DIRS[1:])

    ext = config.add_extension(
        "unuran_wrapper",
        sources=["unuran_wrapper.c"] + sources,
        libraries=[],
        include_dirs=[str(dir_.resolve()) for dir_ in UNURAN_SOURCE_DIRS]
        + [
            os.path.join(
                os.path.dirname(__file__), "..", "..", "_lib", "src"
            )
        ]
        + [os.path.dirname(__file__)],
        language="c",
        define_macros=DEFINE_MACROS,
    )
    ext.pre_build_hook = unuran_pre_build_hook

    config.add_data_files("*.pxd")
    config.add_data_files("*.pyi")

    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration(top_path="").todict())

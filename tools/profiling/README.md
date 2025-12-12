# Flame profiling for SciPy

## Usage

To generate a [flame graph](https://www.brendangregg.com/flamegraphs.html) with https://github.com/nickodell/python-flamegraph,
use:

```console
pixi run flamegraph script_name log_name
```

where `script_name` and `log_name` are such that
the script you want to profile exists at `tools/profiling/script_name.py`,
and `tools/profiling/log_name.log` and `tools/profiling/log_name.svg` will be generated.

To generate a flame graph with py-spy, which can also profile native extensions on Linux, use:

```console
pixi run py-spy script_name log_name
```

> [!CAUTION]
> `py-spy` requires `sudo` on macOS machines.
> You may wish to establish trust in the [the codebase](https://github.com/benfred/py-spy)
> before use.

where `script_name` and `log_name` are such that
the script you want to profile exists at `tools/profiling/script_name.py`,
and `tools/profiling/log_name.log` and `tools/profiling/log_name.svg` will be generated.

## Opportunities for future development

- write this up in the developer documentation
- support for line profiling
- support for Windows?

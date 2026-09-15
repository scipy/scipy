# Flame profiling for SciPy

## Usage

To generate a [flame graph](https://docs.python.org/3.15/library/profiling.sampling.html#flame-graph-format),
use:

```console
pixi run flamegraph script_name output_name
```

where `script_name` and `output_name` are such that
the script you want to profile exists at `tools/profiling/script_name.py`,
and `tools/profiling/output_name.html` will be generated.

## Opportunities for future development

- write this up in the developer documentation
- support for macOS
- support for native frames (with py-spy?)
- support for line profiling
- support for Windows?

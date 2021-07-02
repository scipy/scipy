Enabling SIMD for Pythran
=========================

SIMD support for Pythran can be enabled by using either of the 
following two methods,

* Setting the environment variable, ``export CXXFLAGS='-DUSE_XSIMD'``.
* Adding the argument ``define_macros=[('USE_SIMD', 1)]`` to ``PythranExtension``
  constructor in ``scipy/scipy/interpolate/setup.py``.

Performance Comparison with and without SIMD
--------------------------------------------

The following example, which performs RBF interpolation on images of various sizes,
is used as the basis of comparing computational times with SIMD enabled and disabled::

    import numpy as np
    import pythran, yappi
    import scipy, sys
    from scipy.interpolate import RBFInterpolator
    from tabulate import tabulate
    np.random.seed(0)

    def find_total_and_avg_time(yappi_stats):
        TOTAL_KEY = 6
        AVG_KEY = 14
        total_time = 0
        avg_time = 0
        for stat_dict in yappi_stats:
            total_time += stat_dict[TOTAL_KEY]
            avg_time += stat_dict[AVG_KEY]
        return total_time, avg_time

    print("Pythran Version: ", pythran.__version__)
    print("Scipy Version: ", scipy.__version__)
    print("NumPy Version: ", np.__version__)
    print("Yappi Version: ", '1.3.2')

    results_init = []
    results_eval = []
    flags = '_'.join(sys.argv[1:])
    repeat_freq = 10

    num_points = [40, 50, 60, 70, 80, 90]
    for num in num_points:
        result_init = []
        result_eval = []
        x = np.linspace(0, 1, num)
        y = x[:, None]
        image = x + y
        print(image.shape)

        # Destroy some values
        mask = np.random.random(image.shape) > 0.7
        image[mask] = np.nan

        valid_mask = ~np.isnan(image)
        coords = np.array(np.nonzero(valid_mask)).T
        values = image[valid_mask]
        
        N = coords.shape[0]
        result_init.append(N)
        result_eval.append(N)
        if len(flags) > 0:
            file_name_prefix = 'results/yappi_' + flags + '_'
        else:
            file_name_prefix = 'results/yappi_'
        file = open(file_name_prefix + str(coords.shape[0]) + '_init', 'w')
        with yappi.run(builtins=True):
            for _ in range(repeat_freq):
                it = RBFInterpolator(coords, values)
        init_stats = yappi.get_func_stats()
        init_stats.print_all(out=file)
        total_time, avg_time = find_total_and_avg_time(init_stats)
        result_init.append(total_time)
        result_init.append(avg_time)
        result_init.append(yappi.get_mem_usage())
        file.write("Memory Usage: " + str(yappi.get_mem_usage()) + "\n")
        file.close()
        yappi.clear_stats()

        file = open(file_name_prefix + str(coords.shape[0]) + '_evaluate', 'w')
        with yappi.run(builtins=True):
            for _ in range(repeat_freq):
                filled = it(list(np.ndindex(image.shape))).reshape(image.shape)
        eval_stats = yappi.get_func_stats()
        eval_stats.print_all(out=file)
        total_time, avg_time = find_total_and_avg_time(eval_stats)
        result_eval.append(total_time)
        result_eval.append(avg_time)
        result_eval.append(yappi.get_mem_usage())
        file.write("Memory Usage: " + str(yappi.get_mem_usage()) + "\n")
        file.close()
        yappi.clear_stats()

        del x, y, image, it, filled
        results_init.append(result_init)
        results_eval.append(result_eval)
        print("Completed Successfully for %d x %d image."%(num, num))

    print("Yappi Profiling Summary")
    print("=======================")
    print()
    print("Repeat Frequency:", repeat_freq)
    print()
    print("### Initalising RBFInterpolator\n")
    print(tabulate(results_init, headers=['Number Of Data Points', 
                                        'Total Time (s)', 'Time/Call (s)',
                                        'Total Memory Usage (bytes)'], tablefmt='github'), end="\n\n")
    print("### Evaluating RBFInterpolator\n")
    print(tabulate(results_eval, headers=['Number Of Data Points', 
                                        'Total Time (s)', 'Time/Call (s)',
                                        'Total Memory Usage (bytes)'], tablefmt='github'))

The difference in performance with and without SIMD is negligible as can be seen in the results
below.

Without SIMD
============

Repeat Frequency: 10

### Initalising RBFInterpolator

|   Number Of Data Points |   Total Time (s) |   Time/Call (s) |   Total Memory Usage (bytes) |
|-------------------------|------------------|-----------------|------------------------------|
|                    1102 |          1.20806 |        0.120752 |                       137592 |
|                    1764 |          3.99963 |        0.399886 |                       137592 |
|                    2556 |          7.81915 |        0.781839 |                       137592 |
|                    3453 |         13.9313  |        1.39303  |                       137592 |
|                    4495 |         29.3855  |        2.93842  |                       137592 |
|                    5643 |         50.9823  |        5.09807  |                       137592 |

### Evaluating RBFInterpolator

|   Number Of Data Points |   Total Time (s) |   Time/Call (s) |   Total Memory Usage (bytes) |
|-------------------------|------------------|-----------------|------------------------------|
|                    1102 |          1.085   |        0.102994 |                       137464 |
|                    1764 |          2.49604 |        0.241569 |                       137464 |
|                    2556 |          5.15016 |        0.50252  |                       137464 |
|                    3453 |          9.48217 |        0.931425 |                       137464 |
|                    4495 |         15.9362  |        1.56944  |                       137464 |
|                    5643 |         24.9028  |        2.46105  |                       137464 |

With SIMD
=========

Repeat Frequency: 10

### Initalising RBFInterpolator

|   Number Of Data Points |   Total Time (s) |   Time/Call (s) |   Total Memory Usage (bytes) |
|-------------------------|------------------|-----------------|------------------------------|
|                    1102 |          1.76125 |        0.176073 |                       137592 |
|                    1764 |          4.74023 |        0.473962 |                       137592 |
|                    2556 |          9.40693 |        0.940613 |                       137592 |
|                    3453 |         18.1898  |        1.81887  |                       137592 |
|                    4495 |         33.8257  |        3.38243  |                       137592 |
|                    5643 |         55.6749  |        5.56733  |                       137592 |

### Evaluating RBFInterpolator

|   Number Of Data Points |   Total Time (s) |   Time/Call (s) |   Total Memory Usage (bytes) |
|-------------------------|------------------|-----------------|------------------------------|
|                    1102 |          1.11245 |        0.105558 |                       137464 |
|                    1764 |          2.52389 |        0.244262 |                       137464 |
|                    2556 |          5.61221 |        0.548289 |                       137464 |
|                    3453 |          9.96587 |        0.977694 |                       137464 |
|                    4495 |         17.2801  |        1.70248  |                       137464 |
|                    5643 |         24.7779  |        2.44902  |                       137464 |
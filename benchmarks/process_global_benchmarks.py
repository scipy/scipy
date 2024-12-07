import json
import pandas as pd


def process_global_benchmarks(f):
    """
    Processes the global benchmarks results into pandas DataFrame.

    Parameters
    ----------
    f: {str, file-like}
        Global Benchmarks output

    Returns
    -------
    nfev, success_rate, mean_time
        pd.DataFrame for the mean number of nfev, success_rate, mean_time
        for each optimisation problem.
    """
    with open(f) as fi:
        dct = json.load(fi)

    nfev = []
    nsuccess = []
    mean_time = []

    solvers = dct[list(dct.keys())[0]].keys()
    for problem, results in dct.items():
        _nfev = []
        _nsuccess = []
        _mean_time = []
        for solver, vals in results.items():
            _nfev.append(vals["mean_nfev"])
            _nsuccess.append(vals["nsuccess"] / vals["ntrials"] * 100)
            _mean_time.append(vals["mean_time"])
        nfev.append(_nfev)
        nsuccess.append(_nsuccess)
        mean_time.append(_mean_time)

    nfev = pd.DataFrame(data=nfev, index=dct.keys(), columns=solvers)
    nsuccess = pd.DataFrame(data=nsuccess, index=dct.keys(), columns=solvers)
    mean_time = pd.DataFrame(data=mean_time, index=dct.keys(), columns=solvers)

    return nfev, nsuccess, mean_time

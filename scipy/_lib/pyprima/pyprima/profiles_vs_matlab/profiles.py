import numpy as np
import sys
import os
from excludelist import excludelist
from optiprofiler import set_cutest_problem_options, find_cutest_problems, run_benchmark
import argparse
from time import time
import matlab.engine
import threading
import socket
import struct


nondefault_options = lambda n, f0: {
    'ftarget' : f0 - 314, # if this is doable, 3.14 otherwise
    'maxfev' : 271*n,
    'npt' : int(min(3.14*n, n**1.23)),
    'rhobeg' : 2.718,
    'rhoend' : 3.14*1.0e-4,
}


def get_matlab_options(n, f0):
    if os.environ.get('NONDEFAULT_MATLAB') == 'True':
        options = nondefault_options(n, f0)
        # Change the option name
        options['maxfun'] = matlab.int32(options.pop('maxfev'))
        return options
    return {}


def get_python_options(n, f0):
    if os.environ.get('NONDEFAULT_PYTHON') == 'True':
        return nondefault_options(n, f0)
    return {}


def get_matlab_engine():
    '''
    This function is essentially following the singleton design pattern
    (https://gameprogrammingpatterns.com/singleton.html).
    
    Despite the above link recommending against this pattern we use it for
    two reasons.
    
    1) The "lazy-loading" means that the MATLAB engine won't be started until
    after OptiProfiler has created separate processes via the multiprocessing
    module. If OptiProfiler stops using multiprocessing this approach may need
    to be reviewed.
    
    2) Starting a MATLAB engine takes as much as 30s, so we only want to do it
    once per process and keep the resulting engine. If future versions can
    start up much more quickly then this might not be necessary.
    
    '''
    if not hasattr(get_matlab_engine, 'eng'):
        get_matlab_engine.eng = matlab.engine.start_matlab()
    return get_matlab_engine.eng


def fun_wrapper(server, obj_fun, num_vars):
    server.listen(1)
    conn, _ = server.accept()
    with conn:
        while True:
            bufsize = num_vars * 8  # 8 being sizeof(double)
            flags = socket.MSG_WAITALL  # This helps with large (n ~ 5000) problems, i.e. TESTQUAD
            data = conn.recv(bufsize, flags)
            if not data: break  # This indicates the connection was closed
            if len(data) < bufsize: raise ValueError('Received incomplete data')
            x = struct.unpack(f'{num_vars}d', data)
            x = np.array(x, dtype=np.float64)
            fx = obj_fun(x)
            conn.sendall(struct.pack('d', fx))
            
            
def nl_constraints_wrapper(server, cineq_fun, ceq_fun, num_vars):
    server.listen(1)
    conn, _ = server.accept()
    with conn:
        while True:
            bufsize = num_vars * 8  # 8 being sizeof(double)
            flags = socket.MSG_WAITALL  # This helps with large (n ~ 5000) problems
            data = conn.recv(bufsize, flags)
            if not data: break  # This indicates the connection was closed
            x = struct.unpack(f'{num_vars}d', data)
            x = np.array(x, dtype=np.float64)
            cineq = cineq_fun(x)
            ceq = ceq_fun(x)
            if len(cineq) > 0:
                conn.sendall(struct.pack(f'{len(cineq)}d', *cineq))
            if len(ceq) > 0:
                conn.sendall(struct.pack(f'{len(ceq)}d', *ceq))


def matlab_uobyqa(fun, x0):
    with socket.create_server(('localhost', 0)) as server:
        port = server.getsockname()[1]
        threading.Thread(target=fun_wrapper, args=(server, fun, len(x0))).start()
        f0 = fun.__self__._fun(x0)
        x0_m = matlab.double(x0)
        result = get_matlab_engine().matlab_wrapper(port, 'uobyqa', get_matlab_options(len(x0), f0), x0_m, nargout=4)
        x = np.array(matlab.double(result[0]).tomemoryview().tolist())  # Wrap result in matlab.double in case of n=1
        return x


def python_uobyqa(fun, x0):
    from prima import minimize
    f0 = fun.__self__._fun(x0)
    res = minimize(fun, x0, method='uobyqa', options=get_python_options(len(x0), f0))
    return res.x


def matlab_newuoa(fun, x0):
    # We bind to port 0 so that the OS automatically assigns an available port
    # Then we extract the port number so that we can pass it to the MATLAB engine
    with socket.create_server(('localhost', 0)) as server:
        port = server.getsockname()[1]
        threading.Thread(target=fun_wrapper, args=(server, fun, len(x0))).start()
        # We need to reach into the internals of fun in order to call the original function object
        # so that we don't add to the number of function evaluations. Otherwise OptiProfiler will
        # end up throwing an exception as we hit max_eval. In the future when OptiProfiler offers
        # the signature where it provides the problem itself we'll be able to do this without reaching
        # into the internals like this.
        f0 = fun.__self__._fun(x0)
        x0_m = matlab.double(x0)
        result = get_matlab_engine().matlab_wrapper(port, 'newuoa', get_matlab_options(len(x0), f0), x0_m, nargout=4)
        x = np.array(matlab.double(result[0]).tomemoryview().tolist())  # Wrap result in matlab.double in case of n=1
        return x


def python_newuoa(fun, x0):
    from prima import minimize
    f0 = fun.__self__._fun(x0)
    res = minimize(fun, x0, method='newuoa', options=get_python_options(len(x0), f0))
    return res.x


def matlab_bobyqa(fun, x0, lb, ub):
    with socket.create_server(('localhost', 0)) as server:
        port = server.getsockname()[1]
        threading.Thread(target=fun_wrapper, args=(server, fun, len(x0))).start()
        f0 = fun.__self__._fun(x0)
        x0_m = matlab.double(x0)
        lb = matlab.double(lb)
        ub = matlab.double(ub)
        result = get_matlab_engine().matlab_wrapper(port, 'bobyqa', get_matlab_options(len(x0), f0), x0_m, lb, ub, nargout=4)
        x = np.array(matlab.double(result[0]).tomemoryview().tolist())  # Wrap result in matlab.double in case of n=1
        return x


def python_bobyqa(fun, x0, lb, ub):
    from prima import minimize
    from scipy.optimize import Bounds

    f0 = fun.__self__._fun(x0)
    bounds = Bounds(lb, ub)
    res = minimize(fun, x0, method='bobyqa', bounds=bounds, options=get_python_options(len(x0), f0))
    return res.x


def matlab_lincoa(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq):
    with socket.create_server(('localhost', 0)) as server:
        port = server.getsockname()[1]
        threading.Thread(target=fun_wrapper, args=(server, fun, len(x0))).start()
        f0 = fun.__self__._fun(x0)
        x0_m = matlab.double(x0)
        lb = matlab.double(lb)
        ub = matlab.double(ub)
        a_ub = matlab.double(a_ub)
        b_ub = matlab.double(b_ub)
        a_eq = matlab.double(a_eq)
        b_eq = matlab.double(b_eq)
        result = get_matlab_engine().matlab_wrapper(port, 'lincoa', get_matlab_options(len(x0), f0), x0_m, lb, ub, a_ub, b_ub, a_eq, b_eq, nargout=4)
        x = np.array(matlab.double(result[0]).tomemoryview().tolist())  # Wrap result in matlab.double in case of n=1
        return x


def python_lincoa(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq):
    from prima import minimize, Bounds, LinearConstraint

    f0 = fun.__self__._fun(x0)
    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    res = minimize(fun, x0, method='lincoa', bounds=bounds, constraints=constraints, options=get_python_options(len(x0), f0))
    return res.x


def matlab_cobyla(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ineq, c_eq):
    with socket.create_server(('localhost', 0)) as server:
        port = server.getsockname()[1]
        threading.Thread(target=fun_wrapper, args=(server, fun, len(x0))).start()
        with socket.create_server(('localhost', 0)) as server_nonlcon:
            port_nonlcon = server_nonlcon.getsockname()[1]
            threading.Thread(target=nl_constraints_wrapper, args=(server_nonlcon, c_ineq, c_eq, len(x0))).start()
            f0 = fun.__self__._fun(x0)
            x0_m = matlab.double(x0)
            lb = matlab.double(lb)
            ub = matlab.double(ub)
            a_ub = matlab.double(a_ub)
            b_ub = matlab.double(b_ub)
            a_eq = matlab.double(a_eq)
            b_eq = matlab.double(b_eq)
            c_ineq_x0 = c_ineq(x0)
            m_c_ineq = c_ineq_x0.size
            c_eq_x0 = c_eq(x0)
            m_c_eq = c_eq_x0.size
            result = get_matlab_engine().matlab_wrapper(port, 'cobyla', get_matlab_options(len(x0), f0), x0_m, lb, ub, a_ub, b_ub, a_eq, b_eq, m_c_ineq, m_c_eq, port_nonlcon, nargout=4)
            x = np.array(matlab.double(result[0]).tomemoryview().tolist())  # Wrap result in matlab.double in case of n=1
            return x


def python_cobyla(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ineq, c_eq):
    from prima import minimize, Bounds, LinearConstraint, NonlinearConstraint

    f0 = fun.__self__._fun(x0)
    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    c_ineq_x0 = c_ineq(x0)
    if c_ineq_x0.size > 0:
        constraints.append(NonlinearConstraint(c_ineq, -np.inf, np.zeros_like(c_ineq_x0)))
    c_eq_x0 = c_eq(x0)
    if c_eq_x0.size > 0:
        constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    res = minimize(fun, x0, method='cobyla', bounds=bounds, constraints=constraints, options=get_python_options(len(x0), f0))
    return res.x


def get_problems(description):
    cutest_problem_names = find_cutest_problems(description)
    return list(filter(lambda x: x not in excludelist(description), cutest_problem_names))


if __name__ == '__main__':
    # If we run this script from a directory other than the one that contains it, pycutest's call to importlib will fail,
    # unless we insert the current working directory into the path.
    sys.path.insert(0, os.getcwd())
    os.environ['PYCUTEST_CACHE'] = os.getcwd()
    
    parser = argparse.ArgumentParser(description='Generate performance profiles comparing PRIMA MATLAB to PRIMA Python. '
                                     'By default no profiles are run. To run one or more profiles use the flags '
                                     'specified below. Multiple flags may be specified to run multiple profiles, '
                                     'but please note they will be run sequentially.')
    parser.add_argument('-j', '--n_jobs', type=int, default=None, help='Number of jobs to run in parallel')
    parser.add_argument('-n', '--newuoa', action='store_true', help='Run the NEWUOA benchmark')
    parser.add_argument('-u', '--uobyqa', action='store_true', help='Run the UOBYQA benchmark')
    parser.add_argument('-b', '--bobyqa', action='store_true', help='Run the BOBYQA benchmark')
    parser.add_argument('-l', '--lincoa', action='store_true', help='Run the LINCOA benchmark')
    parser.add_argument('-c', '--cobyla', action='store_true', help='Run the COBYLA benchmark')
    parser.add_argument('--default_only', action='store_true', help='Run only the default options for both MATLAB and Python')
    args = parser.parse_args()
    

    def run_three_benchmarks(matlab_fun, python_fun, algorithm, cutest_problem_names, default_only, n_jobs):
        '''
        Proper validation of both default and nondefault options requires 3 runs: Both default, both nondefault, and
        one default one nondefault. The first two should look identical, and so the third run confirms that our
        experiment setup is valid (i.e. it rules out a scenario where even though options are provided, both algorithms
        end up using default options anyway).
        '''
        # Sharing state with multiprocessing is hard when we can't control the function signature,
        # so we resort to using the environment to pass options.
        os.environ['NONDEFAULT_MATLAB'] = "False"
        os.environ['NONDEFAULT_PYTHON'] = "False"
        algorithm = algorithm.lower()
        ALGORITHM = algorithm.upper()
        project_x0 = algorithm == 'lincoa'
        run_benchmark([matlab_fun, python_fun], [f'MATLAB-{ALGORITHM}', f'Python-{ALGORITHM}'], cutest_problem_names, benchmark_id=f'{algorithm}_default_options', n_jobs=n_jobs, project_x0=project_x0)
        if not default_only:
            os.environ['NONDEFAULT_MATLAB'] = "True"
            os.environ['NONDEFAULT_PYTHON'] = "True"
            run_benchmark([matlab_fun, python_fun], [f'MATLAB-{ALGORITHM}', f'Python-{ALGORITHM}'], cutest_problem_names, benchmark_id=f'{algorithm}_nondefault_options', n_jobs=n_jobs, project_x0=project_x0)
            os.environ['NONDEFAULT_MATLAB'] = "True"
            os.environ['NONDEFAULT_PYTHON'] = "False"
            run_benchmark([matlab_fun, python_fun], [f'MATLAB-{ALGORITHM}', f'Python-{ALGORITHM}'], cutest_problem_names, benchmark_id=f'{algorithm}_different_options', n_jobs=n_jobs, project_x0=project_x0)
    
    if args.newuoa:
        start = time()
        print("Running profiles for NEWUOA")
        with open('uobyqa_newuoa.txt') as f:
            cutest_problem_names = f.read().splitlines()
        cutest_problem_names = list(filter(lambda x: x not in excludelist('unconstrained'), cutest_problem_names))
        run_three_benchmarks(matlab_newuoa, python_newuoa, 'newuoa', cutest_problem_names, args.default_only, args.n_jobs)
        print(f'Completed NEWUOA profile in {time() - start:.2f} seconds')
    
    if args.uobyqa:
        start = time()
        print("Running profiles for UOBYQA")
        with open('uobyqa_newuoa.txt') as f:
            cutest_problem_names = f.read().splitlines()
        cutest_problem_names = list(filter(lambda x: x not in excludelist('unconstrained'), cutest_problem_names))
        run_three_benchmarks(matlab_uobyqa, python_uobyqa, 'uobyqa', cutest_problem_names, args.default_only, args.n_jobs)
        print(f'Completed UOBYQA profile in {time() - start:.2f} seconds')
    
    if args.bobyqa:
        start = time()
        print("Running profiles for BOBYQA")
        with open('bobyqa.txt') as f:
            cutest_problem_names = f.read().splitlines()
        cutest_problem_names = list(filter(lambda x: x not in excludelist('bobyqa'), cutest_problem_names))
        run_three_benchmarks(matlab_bobyqa, python_bobyqa, 'bobyqa', cutest_problem_names, args.default_only, args.n_jobs)
        print(f'Completed BOBYQA profile in {time() - start:.2f} seconds')
    
    if args.lincoa:
        start = time()
        print("Running profiles for LINCOA")
        with open('lincoa.txt') as f:
            cutest_problem_names = f.read().splitlines()
        cutest_problem_names = list(filter(lambda x: x not in excludelist('lincoa'), cutest_problem_names))
        run_three_benchmarks(matlab_lincoa, python_lincoa, 'lincoa', cutest_problem_names, args.default_only, args.n_jobs)
        print(f'Completed LINCOA profile in {time() - start:.2f} seconds')
    
    if args.cobyla:
        start = time()
        print("Running profiles for COBYLA")
        with open('cobyla.txt') as f:
            cutest_problem_names = f.read().splitlines()
        cutest_problem_names = list(filter(lambda x: x not in excludelist('cobyla'), cutest_problem_names))
        run_three_benchmarks(matlab_cobyla, python_cobyla, 'cobyla', cutest_problem_names, args.default_only, args.n_jobs)
        print(f'Completed COBYLA profile in {time() - start:.2f} seconds')

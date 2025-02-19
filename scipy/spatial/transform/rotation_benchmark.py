import os

os.environ["SCIPY_ARRAY_API"] = "1"

import torch
import numpy as np
from typing import Callable
import timeit
import jax
from functools import partial
import jax.numpy as jp
from numpy.typing import NDArray
import matplotlib.pyplot as plt

from typing import Dict


from scipy.spatial.transform import Rotation as R


def create_random_data(
    n_samples: int = 10000, xp_str: str = "numpy", device: str = "cpu"
):
    """Create random test data in numpy format."""
    if xp_str == "numpy":
        return np.random.randn(n_samples, 4), np.random.rand(n_samples, 3)
    elif xp_str == "torch":
        device = "cuda" if device == "gpu" else device
        return torch.rand(n_samples, 4, device=device), torch.rand(
            n_samples, 3, device=device
        )
    elif xp_str == "jax":
        return jax_qp(n_samples, device)
    raise ValueError(f"Invalid xp_str: {xp_str}")


@partial(jax.jit, static_argnums=[0, 1])
def jax_qp(n_samples: int = 10000, device: str = "cpu"):
    dev = jax.devices(device)[0]
    q = jp.array(jax.random.normal(jax.random.PRNGKey(0), (n_samples, 4)), device=dev)
    p = jp.array(jax.random.uniform(jax.random.PRNGKey(0), (n_samples, 3)), device=dev)
    return q, p


def benchmark_function(setup_code: Callable, test_code: Callable) -> NDArray:
    timer = timeit.Timer(stmt=test_code, setup=setup_code)
    R, N = 2, 2
    return np.array(timer.repeat(repeat=R, number=N)) / N


xp_device_combinations = (
    ("numpy", "cpu"),
    ("torch", "cpu"),
    ("torch", "gpu"),
    ("jax", "cpu"),
    ("jax", "gpu"),
)


def benchmark_from_quat(n_samples: int = 10000) -> Dict[str, float]:
    benchmarks = {}

    for xp, device in xp_device_combinations:
        print(f"Benchmarking from_quat with {xp} and {device}")
        q, p, r, from_quat = None, None, None, None

        def setup() -> str:
            nonlocal q, p, r, from_quat
            q, p = create_random_data(n_samples, xp, device)
            dev = "gpu" if "cuda" in str(q.device) else "cpu"
            assert dev == device, f"setup device mismatch: {dev} != {device}"
            if xp == "jax":
                from_quat = jax.jit(partial(R.from_quat, scalar_first=False))
                jax.block_until_ready(from_quat(q))
            r = R.from_quat(q)

        def test():
            nonlocal q
            return R.from_quat(q)

        def jax_test():
            nonlocal q, from_quat
            jax.block_until_ready(from_quat(q))

        timing = benchmark_function(setup, test if xp != "jax" else jax_test)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_as_quat(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark as_quat with different array types."""
    benchmarks = {}

    for xp, device in xp_device_combinations:
        print(f"Benchmarking as_quat with {xp} and {device}")
        q, p, r, as_quat = None, None, None, None

        def setup() -> str:
            nonlocal q, p, r, as_quat
            q, p = create_random_data(n_samples, xp, device)
            dev = "gpu" if "cuda" in str(q.device) else "cpu"
            assert dev == device, f"setup device mismatch: {dev} != {device}"
            r = R.from_quat(q)
            if xp == "jax":
                as_quat = jax.jit(R.as_quat)
                jax.block_until_ready(as_quat(r))

        def test():
            nonlocal r
            return r.as_quat()

        def jax_test():
            nonlocal r, as_quat
            jax.block_until_ready(as_quat(r))

        timing = benchmark_function(setup, test if xp != "jax" else jax_test)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_as_matrix(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark as_matrix with different array types."""
    benchmarks = {}

    for xp, device in xp_device_combinations:
        print(f"Benchmarking as_matrix with {xp} and {device}")
        q, p, r, as_matrix = None, None, None, None

        def setup() -> str:
            nonlocal q, p, r, as_matrix
            q, p = create_random_data(n_samples, xp, device)
            dev = "gpu" if "cuda" in str(q.device) else "cpu"
            assert dev == device, f"setup device mismatch: {dev} != {device}"
            r = R.from_quat(q)
            if xp == "jax":
                as_matrix = jax.jit(R.as_matrix)
                jax.block_until_ready(as_matrix(r))

        def test():
            nonlocal r
            return r.as_matrix()

        def jax_test():
            nonlocal r, as_matrix
            jax.block_until_ready(as_matrix(r))

        timing = benchmark_function(setup, test if xp != "jax" else jax_test)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_apply(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark apply with different array types."""
    benchmarks = {}

    for xp, device in xp_device_combinations:
        print(f"Benchmarking apply with {xp} and {device}")
        q, p, r, apply = None, None, None, None

        def setup() -> str:
            nonlocal q, p, r, apply
            q, p = create_random_data(n_samples, xp, device)
            dev = "gpu" if "cuda" in str(q.device) else "cpu"
            assert dev == device, f"setup device mismatch: {dev} != {device}"
            r = R.from_quat(q)
            if xp == "jax":
                apply = jax.jit(R.apply)
                jax.block_until_ready(apply(r, p))

        def test():
            nonlocal r, p
            return r.apply(p)

        def jax_test():
            nonlocal r, p, apply
            jax.block_until_ready(apply(r, p))

        timing = benchmark_function(setup, test if xp != "jax" else jax_test)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def _benchmark(fn: str, n_samples: int = 10000) -> Dict[str, Dict[str, float]]:
    """Run all rotation benchmarks."""
    if fn == "from_quat":
        return benchmark_from_quat(n_samples)
    elif fn == "as_quat":
        return benchmark_as_quat(n_samples)
    elif fn == "as_matrix":
        return benchmark_as_matrix(n_samples)
    elif fn == "apply":
        return benchmark_apply(n_samples)
    raise ValueError(f"Invalid function: {fn}")


def main():
    sample_sizes = np.logspace(0, 3, 4).astype(int)
    all_results = {}
    fns = ["from_quat", "as_quat", "as_matrix", "apply"]

    for fn in fns:
        all_results[fn] = {}
        for n_samples in sample_sizes:
            print(f"Running {fn} benchmark for {n_samples} samples")
            all_results[fn][n_samples] = _benchmark(fn, n_samples)

    # Plot results

    # Define colors for each XP type and device combination
    colors = {
        "numpy:cpu": "#1f77b4",  # blue
        "torch:cpu": "#ff7f0e",  # orange
        "torch:gpu": "#d62728",  # red
        "jax:cpu": "#2ca02c",  # green
        "jax:gpu": "#9467bd",  # purple
    }

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("Rotation Benchmark Results")
    axes = axes.ravel()

    # Plot one graph per function
    for func_idx, func_name in enumerate(fns):
        ax = axes[func_idx]

        # For each XP type and device combination
        for xp_device in xp_device_combinations:
            xp_device_key = f"{xp_device[0]}:{xp_device[1]}"

            # Collect means and std devs across sample sizes
            means = []
            std_devs = []
            for n_samples in sample_sizes:
                timings = all_results[func_name][n_samples][xp_device_key]
                means.append(np.mean(timings))
                std_devs.append(np.std(timings))

            # Plot with error bars
            ax.errorbar(
                sample_sizes,
                means,
                yerr=std_devs,
                label=xp_device_key,
                color=colors[xp_device_key],
                marker="o",
                capsize=5,
            )

        ax.set_title(func_name)
        ax.set_xlabel("Number of samples")
        ax.set_ylabel("Time (seconds)")
        ax.grid(True)
        ax.legend()

        # Set both axes to logarithmic scale
        ax.set_xscale("log")
        ax.set_yscale("log")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

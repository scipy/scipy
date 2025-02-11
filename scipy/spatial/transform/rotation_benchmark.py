import os

os.environ["SCIPY_ARRAY_API"] = "1"
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"

import torch
import numpy as np
import jax.numpy as jp
from jax.tree_util import register_pytree_node
import timeit
from functools import partial
import jax
from numpy.typing import NDArray
import matplotlib.pyplot as plt

from typing import Dict
from scipy.spatial.transform import Rotation as R

register_pytree_node(R, lambda v: ((v._quat,), None), lambda _, c: R(c[0]))


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
    device = jax.devices(device)[0]
    q = jax.random.normal(jax.random.PRNGKey(0), (n_samples, 4))
    p = jax.random.uniform(jax.random.PRNGKey(0), (n_samples, 3))
    return jax.device_put(q, device), jax.device_put(p, device)


def benchmark_function(setup_code: str, test_code: str) -> NDArray:
    timer = timeit.Timer(stmt=test_code, setup=setup_code)
    R, N = 5, 100
    return np.array(timer.repeat(repeat=R, number=N)) / N


# Common setup code template for benchmarks
SETUP_CODE_TEMPLATE = (
    "import numpy as np\n"
    "import torch\n"
    "import jax.numpy as jp\n"
    "import jax\n"
    "from scipy.spatial.transform import Rotation as R\n"
    "from __main__ import create_random_data\n"
    "q, p = create_random_data({n_samples}, '{xp}', '{device}')\n"
    "r = R.from_quat(q)\n"
    "{extra_setup}"
)
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
        extra_setup = ""
        if xp == "jax":
            extra_setup += "from_quat = jax.jit(R.from_quat)\nfrom_quat(q)\n"
        setup_code = SETUP_CODE_TEMPLATE.format(
            n_samples=n_samples,
            xp=xp,
            device=device,
            extra_setup=extra_setup,
        )
        test_code = (
            "R.from_quat(q)" if xp != "jax" else "jax.block_until_ready(from_quat(q))"
        )
        timing = benchmark_function(setup_code, test_code)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_as_quat(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark as_quat with different array types."""
    benchmarks = {}
    extra_setup = ""

    for xp, device in xp_device_combinations:
        print(f"Benchmarking as_quat with {xp} and {device}")
        if xp == "jax":
            extra_setup += "as_quat = jax.jit(R.as_quat)\nas_quat(r)\n"

        setup_code = SETUP_CODE_TEMPLATE.format(
            n_samples=n_samples,
            xp=xp,
            device=device,
            extra_setup=extra_setup,
        )
        test_code = "r.as_quat()" if xp != "jax" else "as_quat(r).block_until_ready()"
        timing = benchmark_function(setup_code, test_code)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_as_matrix(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark as_matrix with different array types."""
    benchmarks = {}
    extra_setup = ""
    for xp, device in xp_device_combinations:
        print(f"Benchmarking as_matrix with {xp} and {device}")
        if xp == "jax":
            extra_setup += "as_matrix = jax.jit(R.as_matrix)\nas_matrix(r)\n"
        setup_code = SETUP_CODE_TEMPLATE.format(
            n_samples=n_samples,
            xp=xp,
            device=device,
            extra_setup=extra_setup,
        )
        test_code = (
            "r.as_matrix()" if xp != "jax" else "as_matrix(r).block_until_ready()"
        )
        timing = benchmark_function(setup_code, test_code)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_apply(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark apply with different array types."""
    benchmarks = {}
    extra_setup = ""

    for xp, device in xp_device_combinations:
        print(f"Benchmarking apply with {xp} and {device}")
        if xp == "jax":
            extra_setup += "apply = jax.jit(R.apply)\napply(r, p)\n"
        setup_code = SETUP_CODE_TEMPLATE.format(
            n_samples=n_samples,
            xp=xp,
            device=device,
            extra_setup=extra_setup,
        )
        test_code = "r.apply(p)" if xp != "jax" else "apply(r, p).block_until_ready()"
        timing = benchmark_function(setup_code, test_code)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_rotation_functions(n_samples: int = 10000) -> Dict[str, Dict[str, float]]:
    """Run all rotation benchmarks."""
    return {
        "from_quat": benchmark_from_quat(n_samples),
        "as_quat": benchmark_as_quat(n_samples),
        "as_matrix": benchmark_as_matrix(n_samples),
        "apply": benchmark_apply(n_samples),
    }


def main():
    sample_sizes = np.logspace(0, 7, 8).astype(int)
    all_results = {}

    for n_samples in sample_sizes:
        print(f"Running benchmark for {n_samples} samples")
        all_results[n_samples] = benchmark_rotation_functions(n_samples)

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
    for func_idx, func_name in enumerate(
        ["from_quat", "as_quat", "as_matrix", "apply"]
    ):
        ax = axes[func_idx]

        # For each XP type and device combination
        for xp_device in xp_device_combinations:
            xp_device_key = f"{xp_device[0]}:{xp_device[1]}"

            # Collect means and std devs across sample sizes
            means = []
            std_devs = []
            for n_samples in sample_sizes:
                timings = all_results[n_samples][func_name][xp_device_key]
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

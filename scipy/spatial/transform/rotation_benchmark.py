import os

os.environ["SCIPY_ARRAY_API"] = "1"

import torch
import numpy as np
import jax.numpy as jp
from jax.tree_util import register_pytree_node
import timeit
import jax
from numpy.typing import NDArray
import matplotlib.pyplot as plt

from typing import Dict
from scipy.spatial.transform import Rotation as R

register_pytree_node(R, lambda v: ((v._quat,), None), lambda _, c: R(c[0]))


def create_random_data(n_samples: int = 10000):
    """Create random test data in numpy format."""
    np_quat = np.random.rand(n_samples, 4)
    np_quat = np_quat / np.linalg.norm(np_quat, axis=1, keepdims=True)
    np_points = np.random.rand(n_samples, 3)
    return np_quat, np_points


def to_torch(q: np.ndarray, p: np.ndarray, device: str):
    device = "cuda" if device == "gpu" else device
    return torch.tensor(q, device=device), torch.tensor(p, device=device)


def to_jax(q: np.ndarray, p: np.ndarray, device: str):
    device = jax.devices(device)[0]
    return jp.array(q, device=device), jp.array(p, device=device)


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
    "q, p = np.array({quat}), np.array({points})\n"
    "from __main__ import to_torch, to_jax\n"
    "(q, p) = to_torch(q, p, '{device}') if '{xp}' == 'torch' else (q, p)\n"
    "(q, p) = to_jax(q, p, '{device}') if '{xp}' == 'jax' else (q, p)\n"
    "r = R.from_quat(q)\n"
    "{extra_setup}"
)
xp_device_combinations = (
    ("numpy", None),
    ("torch", "cpu"),
    ("torch", "gpu"),
    ("jax", "cpu"),
    ("jax", "gpu"),
)


def benchmark_from_quat(n_samples: int = 10000) -> Dict[str, float]:
    quat, points = create_random_data(n_samples)
    quat = quat.tolist()
    points = points.tolist()
    benchmarks = {}

    for xp, device in xp_device_combinations:
        extra_setup = ""
        if xp == "jax":
            extra_setup += "from_quat = jax.jit(R.from_quat)\nfrom_quat(q)\n"
        setup_code = SETUP_CODE_TEMPLATE.format(
            quat=quat,
            points=points,
            xp=xp,
            device=device,
            extra_setup=extra_setup,
        )
        test_code = (
            "R.from_quat(q)"
            if xp != "jax"
            else "from_quat(q)._quat.block_until_ready()"
        )
        timing = benchmark_function(setup_code, test_code)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_as_quat(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark as_quat with different array types."""
    quat, points = create_random_data(n_samples)
    quat = quat.tolist()
    points = points.tolist()
    benchmarks = {}
    extra_setup = ""

    for xp, device in xp_device_combinations:
        if xp == "jax":
            extra_setup += "as_quat = jax.jit(r.as_quat)\nas_quat()\n"
        setup_code = SETUP_CODE_TEMPLATE.format(
            quat=quat,
            points=points,
            xp=xp,
            device=device,
            extra_setup=extra_setup,
        )
        test_code = "r.as_quat()" if xp != "jax" else "as_quat().block_until_ready()"
        timing = benchmark_function(setup_code, test_code)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_as_matrix(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark as_matrix with different array types."""
    quat, points = create_random_data(n_samples)
    quat = quat.tolist()
    points = points.tolist()
    benchmarks = {}
    extra_setup = ""
    for xp, device in xp_device_combinations:
        if xp == "jax":
            extra_setup += "as_matrix = jax.jit(r.as_matrix)\nas_matrix()\n"
        setup_code = SETUP_CODE_TEMPLATE.format(
            quat=quat,
            points=points,
            xp=xp,
            device=device,
            extra_setup=extra_setup,
        )
        test_code = (
            "r.as_matrix()" if xp != "jax" else "as_matrix().block_until_ready()"
        )
        timing = benchmark_function(setup_code, test_code)
        benchmarks[f"{xp}:{device}"] = timing

    return benchmarks


def benchmark_apply(n_samples: int = 10000) -> Dict[str, float]:
    """Benchmark apply with different array types."""
    quat, points = create_random_data(n_samples)
    quat = quat.tolist()
    points = points.tolist()
    benchmarks = {}
    extra_setup = ""

    for xp, device in xp_device_combinations:
        if xp == "jax":
            extra_setup += "apply = jax.jit(r.apply)\napply(p)\n"
        setup_code = SETUP_CODE_TEMPLATE.format(
            quat=quat,
            points=points,
            xp=xp,
            device=device,
            extra_setup=extra_setup,
        )
        test_code = "r.apply(p)" if xp != "jax" else "apply(p).block_until_ready()"
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
    sample_sizes = np.logspace(0, 4, 5).astype(int)
    all_results = {}

    for n_samples in sample_sizes:
        all_results[n_samples] = benchmark_rotation_functions(n_samples)

    # Plot results

    # Define colors for each XP type and device combination
    colors = {
        "numpy:None": "#1f77b4",  # blue
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

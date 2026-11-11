#include "HConfig.h"
#include "gpu_backend.hpp"
#include "pdhg_kernels.hpp"
#include <cmath>
#include <cstdio>

#ifdef HIGHSINT64
typedef int64_t HighsInt;
typedef uint64_t HighsUInt;
#define HIGHSINT_FORMAT PRId64
#else
typedef int HighsInt;
typedef unsigned int HighsUInt;
#define HIGHSINT_FORMAT "d"
#endif

// Define Infinity for GPU
#define GPU_INF 1e20

// Buffer Indices for the reduction array
#define IDX_PRIMAL_FEAS 0
#define IDX_DUAL_FEAS   1
#define IDX_PRIMAL_OBJ  2
#define IDX_DUAL_OBJ    3

// warpSize is a built-in device variable that returns the actual hardware
// warp/wavefront width (32 on CUDA and RDNA, 64 on GCN/CDNA).
__device__ __forceinline__ double gpuReduce(double val) {
  for (int offset = warpSize / 2; offset > 0; offset >>= 1) {
#ifdef USE_HIP
    val += __shfl_down(val, offset);
#else
    val += __shfl_down_sync(0xffffffff, val, offset);
#endif
  }
  return val;
}

// Utility for robust 1D kernel launches
#define GPU_GRID_STRIDE_LOOP(i, n)                                  \
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < (n); \
       i += blockDim.x * gridDim.x)

// Helper function to calculate launch configuration
static dim3 GetLaunchConfig(int n, int block_size = 256) {
  int num_blocks = (n + block_size - 1) / block_size;
  return dim3(num_blocks, 1, 1);
}

// ============================================================================
// FMA (Fused Multiply-Add) Helper
// ============================================================================
// Use __fma_rn for different precision:
//   __fma_rn  - double precision, round to nearest
//   __fmaf_rn - single precision, round to nearest
// FMA computes (a * b + c) in a single operation with only one rounding,
// which is faster and more numerically accurate.

__device__ __forceinline__ double fma_rn(double a, double b, double c) {
    return __fma_rn(a, b, c);  // a * b + c
}

// === KERNEL 1: Update X (Primal Step) ===
__global__ void kernelUpdateX(
    double* __restrict__ d_x_new, const double* __restrict__ d_x_old,
    const double* __restrict__ d_aty, const double* __restrict__ d_cost,
    const double* __restrict__ d_lower, const double* __restrict__ d_upper,
    double primal_step, int n_cols)
{
  GPU_GRID_STRIDE_LOOP(i, n_cols) {
    double x_updated = fma_rn(primal_step, d_aty[i] - d_cost[i], d_x_old[i]);

    // Project to bounds [l, u]
    d_x_new[i] = fmin(fmax(x_updated, d_lower[i]), d_upper[i]);
  }
}

// === KERNEL 2: Update Y (Dual Step) ===
__global__ void kernelUpdateY(
    double* __restrict__ d_y_new, const double* __restrict__ d_y_old,
    const double* __restrict__ d_ax_old, const double* __restrict__ d_ax_new,
    const double* __restrict__ d_rhs, const bool* __restrict__ d_is_equality,
    double dual_step, int n_rows)
{
  GPU_GRID_STRIDE_LOOP(j, n_rows) {
    double residual = fma_rn(-2.0, d_ax_new[j], d_rhs[j] + d_ax_old[j]);
    double dual_update = fma_rn(dual_step, residual, d_y_old[j]);
    d_y_new[j] = d_is_equality[j] ? dual_update : fmax(0.0, dual_update);
  }
}

// === KERNEL 3: Update Averages ===
// x_sum = x_sum + weight * x_next.
// y_sum = y_sum + weight * y_next
__global__ void kernelUpdateAverages(
    double* __restrict__ d_x_sum, double* __restrict__ d_y_sum,
    const double* __restrict__ d_x_next, const double* __restrict__ d_y_next,
    double weight, int n_cols, int n_rows)
{
    // Update x_sum
    GPU_GRID_STRIDE_LOOP(i, n_cols) {
        d_x_sum[i] = fma_rn(weight, d_x_next[i], d_x_sum[i]);
    }

    // Update y_sum
    GPU_GRID_STRIDE_LOOP(j, n_rows) {
        d_y_sum[j] = fma_rn(weight, d_y_next[j], d_y_sum[j]);
    }
}

__global__ void kernelScaleVector(
    double* __restrict__ d_out, const double* __restrict__ d_in,
    double scale, int n)
{
  GPU_GRID_STRIDE_LOOP(i, n) {
    d_out[i] = d_in[i] * scale;
  }
}

// === KERNEL 4: Primal Convergence Check (Row-wise) ===
__global__ void kernelCheckPrimal(
  double* __restrict__ d_results,
  const double* __restrict__ d_ax, const double* __restrict__ d_y,
  const double* __restrict__ d_row_lower, const double* d_row_scale,
  const bool* __restrict__ d_is_equality, int n_rows){
  double local_feas_sq = 0.0;
  double local_dual_obj = 0.0;

  GPU_GRID_STRIDE_LOOP(i,n_rows){
    double val_y = d_y[i];
    double val_b = d_row_lower[i];
    double val_ax = d_ax[i];

    if (abs(val_b) < GPU_INF){
      local_dual_obj = fma_rn(val_b, val_y, local_dual_obj);
    }

    double residual = val_ax - val_b;

    if (!d_is_equality[i]) {
      residual = fmin(0.0, residual);
    }

    if (d_row_scale != nullptr) {
      residual *= d_row_scale[i];
    }

    local_feas_sq = fma_rn(residual, residual, local_feas_sq);
  }

  // Warp-level reduction
  local_feas_sq  = gpuReduce(local_feas_sq);
  local_dual_obj = gpuReduce(local_dual_obj);

  // Block-level reduction: one atomicAdd per block to global memory instead of
  // one per warp, cutting contention by ~blockDim.x/warpSize.
  __shared__ double s_feas[2];
  if (threadIdx.x == 0) { s_feas[0] = 0.0; s_feas[1] = 0.0; }
  __syncthreads();

  if ((threadIdx.x & (warpSize - 1)) == 0) {
    atomicAdd(&s_feas[0], local_feas_sq);
    atomicAdd(&s_feas[1], local_dual_obj);
  }
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(&d_results[IDX_PRIMAL_FEAS], s_feas[0]);
    atomicAdd(&d_results[IDX_DUAL_OBJ],    s_feas[1]);
  }
}

// === KERNEL 5: Dual Convergence Check (Column-wise) ===
// d_slack_pos doubles as input (raw Halpern slack) and output (s_pos).
__global__ void kernelCheckDual(
    double* __restrict__ d_results,
    double* __restrict__ d_slack_pos,   // Input (raw halpern slack) / Output (s_pos)
    double* __restrict__ d_slack_neg,   // Output
    const double* __restrict__ d_aty,
    const double* __restrict__ d_x,
    const double* __restrict__ d_cost,      // c
    const double* __restrict__ d_col_lower, // l
    const double* __restrict__ d_col_upper, // u
    const double* d_col_scale,              // Can be nullptr
    int n_cols,
    bool use_halpern_slack)
{
  double local_dual_feas_sq = 0.0;
  double local_primal_obj = 0.0;
  double local_dual_obj_part = 0.0;

  GPU_GRID_STRIDE_LOOP(i, n_cols){
    double val_x = d_x[i];
    double val_c = d_cost[i];
    double val_aty = d_aty[i];
    double val_l = d_col_lower[i];
    double val_u = d_col_upper[i];

    local_primal_obj = fma_rn(val_c, val_x, local_primal_obj);
    double dual_residual = val_c - val_aty;

    double s_pos = 0.0;
    double s_neg = 0.0;

    if (use_halpern_slack) {
      // Read the raw Halpern slack previously saved in d_slack_pos by primal major step
      double raw_slack = d_slack_pos[i];
      s_pos = fmax(0.0, raw_slack);
      s_neg = fmax(0.0, -raw_slack);
    } else {
      // Standard bounds-based projection (for Average iterates)
      if (val_l > -GPU_INF) s_pos = fmax(0.0, dual_residual);
      if (val_u < GPU_INF)  s_neg = fmax(0.0, -dual_residual);
    }

    // Write back the clean positive/negative parts for objective calc & host copy
    d_slack_pos[i] = s_pos;
    d_slack_neg[i] = s_neg;

    double eff_dual_residual = dual_residual - s_pos + s_neg;

    if (d_col_scale != nullptr) {
      eff_dual_residual *= d_col_scale[i];
    }

    local_dual_feas_sq = fma_rn(eff_dual_residual, eff_dual_residual, local_dual_feas_sq);

    if (val_l > -GPU_INF) local_dual_obj_part = fma_rn(val_l,  s_pos, local_dual_obj_part);
    if (val_u <  GPU_INF) local_dual_obj_part = fma_rn(-val_u, s_neg, local_dual_obj_part);
  }

  // Warp-level reduction
  local_dual_feas_sq  = gpuReduce(local_dual_feas_sq);
  local_primal_obj    = gpuReduce(local_primal_obj);
  local_dual_obj_part = gpuReduce(local_dual_obj_part);

  // Block-level reduction via shared memory before global atomicAdd.
  __shared__ double s_dual[3];
  if (threadIdx.x == 0) { s_dual[0] = 0.0; s_dual[1] = 0.0; s_dual[2] = 0.0; }
  __syncthreads();

  if ((threadIdx.x & (warpSize - 1)) == 0) {
    atomicAdd(&s_dual[0], local_dual_feas_sq);
    atomicAdd(&s_dual[1], local_primal_obj);
    atomicAdd(&s_dual[2], local_dual_obj_part);
  }
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(&d_results[IDX_DUAL_FEAS],  s_dual[0]);
    atomicAdd(&d_results[IDX_PRIMAL_OBJ], s_dual[1]);
    atomicAdd(&d_results[IDX_DUAL_OBJ],   s_dual[2]);
  }
}

// ============================================================================
// HALPERN-MODE
// ============================================================================

// Minor: every iteration. Only writes reflected_primal.
__global__ void kernelHalpernPrimalMinor(
    const double* __restrict__ current_primal,
    double* __restrict__ reflected_primal,
    const double* __restrict__ dual_product,     // A^T * current_dual
    const double* __restrict__ objective,
        const double* __restrict__ var_lb,
        const double* __restrict__ var_ub,
        const double* __restrict__ step_size_ptr, int n)
{
        const double step_size = *step_size_ptr;
    GPU_GRID_STRIDE_LOOP(i, n) {
        double temp = fma_rn(-step_size, objective[i] - dual_product[i],
                             current_primal[i]);
        double temp_proj = fmax(var_lb[i], fmin(temp, var_ub[i]));
        reflected_primal[i] = fma_rn(2.0, temp_proj, -current_primal[i]);
    }
}

// Major: convergence-check iterations.
// Additionally writes pdhg_primal and dual_slack.
__global__ void kernelHalpernPrimalMajor(
    const double* __restrict__ current_primal,
    double* __restrict__ pdhg_primal,
    double* __restrict__ reflected_primal,
    const double* __restrict__ dual_product,
    const double* __restrict__ objective,
        const double* __restrict__ var_lb,
        const double* __restrict__ var_ub,
        const double* __restrict__ step_size_ptr, int n,
    double* __restrict__ dual_slack)
{
        const double step_size = *step_size_ptr;
        const double inv_step_size = 1.0 / step_size;
    GPU_GRID_STRIDE_LOOP(i, n) {
        double temp = fma_rn(-step_size, objective[i] - dual_product[i],
                             current_primal[i]);
        double temp_proj = fmax(var_lb[i], fmin(temp, var_ub[i]));
        pdhg_primal[i] = temp_proj;
        dual_slack[i] = (temp_proj - temp) * inv_step_size;
        reflected_primal[i] = fma_rn(2.0, temp_proj, -current_primal[i]);
    }
}

__global__ void kernelHalpernDualMinor(
    const double* __restrict__ current_dual,
    double* __restrict__ reflected_dual,
    const double* __restrict__ primal_product,   // A * reflected_primal
    const double* __restrict__ rhs,              // row_lower (b)
    const bool* __restrict__ is_equality,
        const double* __restrict__ step_size_ptr, int n)
{
        const double step_size = *step_size_ptr;
    GPU_GRID_STRIDE_LOOP(i, n) {
        double dual_update = fma_rn(step_size, rhs[i] - primal_product[i],
                                    current_dual[i]);
        double pdhg_dual_i = is_equality[i] ? dual_update
                                             : fmax(0.0, dual_update);
        reflected_dual[i] = fma_rn(2.0, pdhg_dual_i, -current_dual[i]);
    }
}

__global__ void kernelHalpernDualMajor(
    const double* __restrict__ current_dual,
    double* __restrict__ pdhg_dual,
    double* __restrict__ reflected_dual,
    const double* __restrict__ primal_product,
    const double* __restrict__ rhs,
    const bool* __restrict__ is_equality,
        const double* __restrict__ step_size_ptr, int n)
{
        const double step_size = *step_size_ptr;
    GPU_GRID_STRIDE_LOOP(i, n) {
        double dual_update = fma_rn(step_size, rhs[i] - primal_product[i],
                                    current_dual[i]);
        double pdhg_dual_i = is_equality[i] ? dual_update
                                             : fmax(0.0, dual_update);
        pdhg_dual[i] = pdhg_dual_i;
        reflected_dual[i] = fma_rn(2.0, pdhg_dual_i, -current_dual[i]);
    }
}

// current = w * (ρ*reflected + (1-ρ)*current) + (1-w)*initial
__global__ void kernelHalpernBlend(
    double* __restrict__ current,
    const double* __restrict__ reflected,
    const double* __restrict__ initial,
    const int* __restrict__ halpern_iteration,
    int k_offset,
    double reflection_coeff,   // ρ, typically 1.0
    HighsInt n)
{
    const int current_k           = (*halpern_iteration) + k_offset;
    const double weight           = static_cast<double>(current_k) / (current_k + 1.0);
    const double one_minus_weight = 1.0 - weight;
    const double one_minus_rho    = 1.0 - reflection_coeff;
    GPU_GRID_STRIDE_LOOP(i, n) {
        double blended = fma_rn(reflection_coeff, reflected[i],
                                fma_rn(one_minus_rho, current[i], 0.0));
        current[i]     = fma_rn(weight, blended,
                                fma_rn(one_minus_weight, initial[i], 0.0));
    }
}

// Add C++ wrapper functions to launch the kernels
extern "C" {
void launchKernelUpdateX_wrapper(
    double* d_x_new, const double* d_x_old, const double* d_aty,
    const double* d_cost, const double* d_lower, const double* d_upper,
    double primal_step, int n_cols, gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config = GetLaunchConfig(n_cols, block_size);

    kernelUpdateX<<<config.x, block_size, 0, stream>>>(
        d_x_new, d_x_old, d_aty,
        d_cost, d_lower, d_upper,
        primal_step, n_cols);

    gpuGetLastError();
}

void launchKernelUpdateY_wrapper(
    double* d_y_new, const double* d_y_old,
    const double* d_ax_old, const double* d_ax_new,
    const double* d_rhs, const bool* d_is_equality,
    double dual_step, int n_rows, gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config = GetLaunchConfig(n_rows, block_size);

    kernelUpdateY<<<config.x, block_size, 0, stream>>>(
        d_y_new, d_y_old,
        d_ax_old, d_ax_new,
        d_rhs, d_is_equality,
        dual_step, n_rows);

    gpuGetLastError();
}

void launchKernelUpdateAverages_wrapper(
    double* d_x_sum, double* d_y_sum,
    const double* d_x_next, const double* d_y_next,
    double weight, int n_cols, int n_rows, gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config_x = GetLaunchConfig(n_cols, block_size);
    dim3 config_y = GetLaunchConfig(n_rows, block_size);
    kernelUpdateAverages<<<config_x.x > config_y.x ? config_x.x : config_y.x, block_size, 0, stream>>>(
        d_x_sum, d_y_sum,
        d_x_next, d_y_next,
        weight, n_cols, n_rows);
    gpuGetLastError();
}

void launchKernelScaleVector_wrapper(
    double* d_out, const double* d_in,
    double scale, int n, gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config = GetLaunchConfig(n, block_size);

    kernelScaleVector<<<config.x, block_size, 0, stream>>>(
        d_out, d_in, scale, n);

    gpuGetLastError();
}

void launchCheckConvergenceKernels_wrapper(
    double* d_results,
    double* d_slack_pos, double* d_slack_neg,
    const double* d_x, const double* d_y,
    const double* d_ax, const double* d_aty,
    const double* d_col_cost, const double* d_row_lower,
    const double* d_col_lower, const double* d_col_upper,
    const bool* d_is_equality,
    const double* d_col_scale, const double* d_row_scale,
    int n_cols, int n_rows,
    bool use_halpern_slack,
    gpuStream_t stream)
{
    // 1. Zero out results
    gpuMemsetAsync(d_results, 0, 4 * sizeof(double), stream);

    int block_size = 256;

    // 2. Launch Primal Kernel
    dim3 grid_rows = GetLaunchConfig(n_rows, block_size);
    kernelCheckPrimal<<<grid_rows, block_size, 0, stream>>>(
        d_results, d_ax, d_y, d_row_lower, d_row_scale, d_is_equality, n_rows
    );

    // 3. Launch Dual Kernel
    dim3 grid_cols = GetLaunchConfig(n_cols, block_size);
    kernelCheckDual<<<grid_cols, block_size, 0, stream>>>(
        d_results, d_slack_pos, d_slack_neg, d_aty, d_x,
        d_col_cost, d_col_lower, d_col_upper, d_col_scale, n_cols,
        use_halpern_slack
    );

    gpuGetLastError();
}

void launchKernelHalpernPrimalMinor_wrapper(
    const double* d_current_primal, double* d_reflected_primal,
    const double* d_dual_product, const double* d_objective,
    const double* d_var_lb, const double* d_var_ub,
    const double* d_step_size, int n, gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config = GetLaunchConfig(n, block_size);
    kernelHalpernPrimalMinor<<<config.x, block_size, 0, stream>>>(
        d_current_primal, d_reflected_primal, d_dual_product,
        d_objective, d_var_lb, d_var_ub, d_step_size, n);
    gpuGetLastError();
}

void launchKernelHalpernPrimalMajor_wrapper(
    const double* d_current_primal, double* d_pdhg_primal,
    double* d_reflected_primal, const double* d_dual_product,
    const double* d_objective, const double* d_var_lb,
    const double* d_var_ub, const double* d_step_size, int n,
    double* d_dual_slack, gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config = GetLaunchConfig(n, block_size);
    kernelHalpernPrimalMajor<<<config.x, block_size, 0, stream>>>(
        d_current_primal, d_pdhg_primal, d_reflected_primal,
        d_dual_product, d_objective, d_var_lb, d_var_ub,
        d_step_size, n, d_dual_slack);
    gpuGetLastError();
}

void launchKernelHalpernDualMinor_wrapper(
    const double* d_current_dual, double* d_reflected_dual,
    const double* d_primal_product, const double* d_rhs,
    const bool* d_is_equality, const double* d_step_size, int n,
    gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config = GetLaunchConfig(n, block_size);
    kernelHalpernDualMinor<<<config.x, block_size, 0, stream>>>(
        d_current_dual, d_reflected_dual, d_primal_product,
        d_rhs, d_is_equality, d_step_size, n);
    gpuGetLastError();
}

void launchKernelHalpernDualMajor_wrapper(
    const double* d_current_dual, double* d_pdhg_dual,
    double* d_reflected_dual, const double* d_primal_product,
    const double* d_rhs, const bool* d_is_equality,
    const double* d_step_size, int n, gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config = GetLaunchConfig(n, block_size);
    kernelHalpernDualMajor<<<config.x, block_size, 0, stream>>>(
        d_current_dual, d_pdhg_dual, d_reflected_dual,
        d_primal_product, d_rhs, d_is_equality, d_step_size, n);
    gpuGetLastError();
}

void launchKernelHalpernBlend_wrapper(
    double* d_current, const double* d_reflected,
    const double* d_initial, const int* d_halpern_iteration, int k_offset,
    double reflection_coeff, int n, gpuStream_t stream)
{
    const int block_size = 256;
    dim3 config = GetLaunchConfig(n, block_size);
    kernelHalpernBlend<<<config.x, block_size, 0, stream>>>(
        d_current, d_reflected, d_initial,
        d_halpern_iteration, k_offset, reflection_coeff, n);
    gpuGetLastError();
}
} // extern "C"

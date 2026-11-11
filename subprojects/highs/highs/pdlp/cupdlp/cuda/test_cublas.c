#include "cupdlp_cuda_kernels.cuh"
#include "cupdlp_cudalinalg.cuh"

int use_cublas(cublasHandle_t cublashandle) {
  cupdlp_int len = 10;
  // cupdlp_int len = 1<<10;

  // int N = 1<<20;

  // alloc and init host vec memory
  cupdlp_float *h_vec1 = (cupdlp_float *)malloc(len * sizeof(cupdlp_float));
  cupdlp_float *h_vec2 = (cupdlp_float *)malloc(len * sizeof(cupdlp_float));
  for (cupdlp_int i = 0; i < len; i++) {
    h_vec1[i] = 1.0;
    h_vec2[i] = i;
    // h_vec1[i] = 1.0;
    // h_vec2[i] = 2.0;
  }

  // alloc and init device vec memory
  cupdlp_float *d_vec1;
  cupdlp_float *d_vec2;
  CHECK_CUDA(cudaMalloc((void **)&d_vec1, len * sizeof(cupdlp_float)))
  // CHECK_CUDA(cudaMemcpy(d_vec1, h_vec1, len * sizeof(cupdlp_float),
  //                       cudaMemcpyHostToDevice))

  CHECK_CUDA(cudaMalloc((void **)&d_vec2, len * sizeof(cupdlp_float)))
  CHECK_CUDA(cudaMemcpy(d_vec2, h_vec2, len * sizeof(cupdlp_float),
                        cudaMemcpyHostToDevice))

  // init cublas handle
  // cublasHandle_t cublashandle;
  // CHECK_CUBLAS(cublasCreate(&cublashandle));

  cupdlp_float result;
  // call nrm2
  CHECK_CUBLAS(cublasDnrm2(cublashandle, len, d_vec1, 1, &result));

  // print result
  printf("2-norm is :%f\n", result);

  // copy result back to host
  // cudaMemcpy(h_vec1, d_vec1, len * sizeof(cupdlp_float),
  //            cudaMemcpyDeviceToHost);
  // cudaMemcpy(h_vec2, d_vec2, len * sizeof(cupdlp_float),
  //            cudaMemcpyDeviceToHost);
  // cudaError_t errSync = cudaGetLastError();
  // cudaError_t errAsync = cudaDeviceSynchronize();
  // if (errSync != cudaSuccess)
  //     printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
  // if (errAsync != cudaSuccess)
  //     printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));

  // // print result
  // for (cupdlp_int i = 0; i < len; i++) {
  //     printf("%f\n", h_vec1[i]);
  //     // printf("%f\n", h_vec2[i]);
  // }

  // destroy cublas handle
  //  CHECK_CUBLAS(cublasDestroy(cublashandle));

  // free memory
  free(h_vec1);
  free(h_vec2);
  CHECK_CUDA(cudaFree(d_vec1))
  CHECK_CUDA(cudaFree(d_vec2))

  return 0;
}
int main() {
  // try cupdlp_edot_cuda

  // int nDevices;

  // cudaGetDeviceCount(&nDevices);
  //    for (int i = 0; i < nDevices; i++) {
  //        cudaDeviceProp prop;
  //        cudaGetDeviceProperties(&prop, i);
  //        printf("Device Number: %d\n", i);
  //        printf("  Device name: %s\n", prop.name);
  //        printf("  Memory Clock Rate (KHz): %d\n",
  //               prop.memoryClockRate);
  //        printf("  Memory Bus Width (bits): %d\n",
  //               prop.memoryBusWidth);
  //        printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
  //               2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8)
  //               / 1.0e6);
  //    }

  // cupdlp_int len = 10;
  // // cupdlp_int len = 1<<10;

  // // int N = 1<<20;

  // // alloc and init host vec memory
  // cupdlp_float *h_vec1 = (cupdlp_float *) malloc(len * sizeof(cupdlp_float));
  // cupdlp_float *h_vec2 = (cupdlp_float *) malloc(len * sizeof(cupdlp_float));
  // for (cupdlp_int i = 0; i < len; i++) {
  //     h_vec1[i] = 1.0;
  //     h_vec2[i] = i;
  //     // h_vec1[i] = 1.0;
  //     // h_vec2[i] = 2.0;
  // }

  // // alloc and init device vec memory
  // cupdlp_float *d_vec1;
  // cupdlp_float *d_vec2;
  // cudaMalloc((void **) &d_vec1, len * sizeof(cupdlp_float));
  // // cudaMemcpy(d_vec1, h_vec1, len * sizeof(cupdlp_float),
  // //            cudaMemcpyHostToDevice);

  // cudaMalloc((void **) &d_vec2, len * sizeof(cupdlp_float));
  // cudaMemcpy(d_vec2, h_vec2, len * sizeof(cupdlp_float),
  //            cudaMemcpyHostToDevice);

  // init cublas handle
  cublasHandle_t cublashandle;
  CHECK_CUBLAS(cublasCreate(&cublashandle));
  use_cublas(cublashandle);
  // cupdlp_float result;
  // call nrm2
  // CHECK_CUBLAS(cublasDnrm2(cublashandle, len, d_vec1, 1, &result));

  // print result
  // printf("2-norm is :%f\n", result);

  // copy result back to host
  // cudaMemcpy(h_vec1, d_vec1, len * sizeof(cupdlp_float),
  //            cudaMemcpyDeviceToHost);
  // cudaMemcpy(h_vec2, d_vec2, len * sizeof(cupdlp_float),
  //            cudaMemcpyDeviceToHost);
  // cudaError_t errSync = cudaGetLastError();
  // cudaError_t errAsync = cudaDeviceSynchronize();
  // if (errSync != cudaSuccess)
  //     printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
  // if (errAsync != cudaSuccess)
  //     printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));

  // // print result
  // for (cupdlp_int i = 0; i < len; i++) {
  //     printf("%f\n", h_vec1[i]);
  //     // printf("%f\n", h_vec2[i]);
  // }

  // destroy cublas handle
  CHECK_CUBLAS(cublasDestroy(cublashandle));

  // // free memory
  // free(h_vec1);
  // free(h_vec2);
  // cudaFree(d_vec1);
  // cudaFree(d_vec2);
}
#ifndef PROPACK__COMMON_H
#define PROPACK__COMMON_H

void scompute_mu(float* restrict mu, const int j, const float delta, const float eta, int* restrict indices);
void sset_mu(const int k, float* restrict mu, int* const restrict indices, const float val);
void supdate_mu(float* mumax, float* restrict mu, float* restrict nu, const int j, float* restrict alpha, float* restrict beta, const float anorm, const float eps1);
void supdate_nu(float* numax, float* restrict mu, float* restrict nu, const int j, float* restrict alpha, float* restrict beta, const float anorm, const float eps1);

void dcompute_mu(double* restrict mu, const int j, const double delta, const double eta, int* restrict indices);
void dset_mu(const int k, double* restrict mu, int* const restrict indices, const double val);
void dupdate_mu(double* mumax, double* restrict mu, double* restrict nu, const int j, double* restrict alpha, double* restrict beta, const double anorm, const double eps1);
void dupdate_nu(double* numax, double* restrict mu, double* restrict nu, const int j, double* restrict alpha, double* restrict beta, const double anorm, const double eps1);


#endif // PROPACK__COMMON_H

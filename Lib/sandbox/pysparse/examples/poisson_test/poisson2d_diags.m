function L = poisson2d_diags(n)
  e = ones(n*n,1);
  f = e;
  f(n:n:n*n) = 0;
  L = spdiags([-e -f 4*e -f -e], [-n -1 0 1 n], n*n, n*n);

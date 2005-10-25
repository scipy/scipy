function L = poisson2d_kron(n)
  e = ones(n,1);
  P = spdiags([-e 2*e -e], [-1 0 1], n, n);
  L = kron(P,speye(n)) + kron(speye(n),P);

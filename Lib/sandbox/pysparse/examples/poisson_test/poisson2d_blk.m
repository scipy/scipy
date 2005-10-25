function L = poisson2d_blk(n)
  e = ones(n,1);
  P = spdiags([-e 4*e -e], [-1 0 1], n, n);
  I = -speye(n);
  L = sparse(n*n);
  for i = 1:n:n*n
    L(i:i+n-1,i:i+n-1) = P;
    if i > 1, L(i:i+n-1,i-n:i-1) = I; end
    if i < n*n - n, L(i:i+n-1,i+n:i+2*n-1) = I; end
  end

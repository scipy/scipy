% Create sparse 2D Poisson matrix
%
% Used for benchmarking element-wise assignment to sparse matrices in Matlab.
%
% $Id: poisson2d.m,v 1.3 2003/11/17 13:22:58 geus Exp $

function L = poisson2d(n)
L = spalloc(n*n, n*n, 5*n*n);
  for i = 1:n
    for j = 1:n
      k = i + n*(j-1);
      L(k,k) = 4;
      if i > 1, L(k,k-1) = -1; end
      if i < n, L(k,k+1) = -1; end
      if j > 1, L(k,k-n) = -1; end
      if j < n, L(k,k+n) = -1; end
    end
  end

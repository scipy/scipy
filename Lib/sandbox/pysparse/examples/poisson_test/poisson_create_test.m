% Benchmark for measuring performance of creating sparse matrices 
% using single element assignments.
%
% Used for comparison with Python and native C/C++
%
%
% Results for Matlab 6.5 (R13) on nedelec (Intel P4 2.40GHz, 512 KB cache)
%
%        n        t_constr
%    ---------------------
%      100          4.5055
%      300        608.5301
%      500       4648.1
%     1000             inf
%
% $Id: poisson_create_test.m,v 1.1 2003/11/17 15:40:14 geus Exp $

tic
A = poisson2d(100);
nnz(A)
toc

tic
A = poisson2d(300);
nnz(A)
toc

tic
A = poisson2d(500);
nnz(A)
toc

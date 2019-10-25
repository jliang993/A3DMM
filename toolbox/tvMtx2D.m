function [Dx, Dy] = tvMtx2D(n)
% For SQUARE mtx first
% should be sparse

A = -diag(ones(n,1))+diag(ones(n-1,1),1);
% A(end,:) = 0;

E = eye(n);

Dx = kron(A, E);
Dy = kron(E, A);
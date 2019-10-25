function x = vecsvt(mtx, tau, n)
% input mtx is a column vector 
% for economy svd, s is square

% n = sqrt(max(size(mtx)));

mtx = reshape(mtx, n);

[u,s,v] = svd(mtx, 'econ');

s = diag( wthresh( diag(s), 's', tau ) );

x = u*s*(v');
x = x(:);


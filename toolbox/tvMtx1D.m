function [D] = tvMtx1D(n)

D = -diag(ones(n,1))+diag(ones(n-1,1),1);
D(end,end) = 0;


% D = diag(ones(n,1)) - diag(ones(n-1,1),1);
% D(1,1) = 0;
% 
% D = D';
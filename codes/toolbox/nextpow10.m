function N = nextpow10(n)
%NEXTPOW10 Next higher power of 10.
%   NEXTPOW10(N) returns the first P such that 10.^P >= abs(N).
%

p = log10(abs(n));

% f1 = floor(p);
% f2 = p - f1;
% N = ceil(10^f2)*10^f1;

f1 = ceil(p);
N = 10^(f1);
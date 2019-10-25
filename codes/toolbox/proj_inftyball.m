function y = proj_inftyball(x, gamma)

if nargin<2
    gamma = 1;
end

y = x;

idx = abs(x) >= gamma;

y(idx) = gamma* sign(x(idx));
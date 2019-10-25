function g = grad_logistic(x, Wy)

% [m, n] = size(W);


% Wy = Y .* W;

e = exp(-Wy*x);
w = ( 1 ./ (1 + e) - 1 );
% Wn = repmat(W, 1, n);

g = (Wy')* w;

% P = Wn .* K;
% g = ( sum(P) - sum(K) )';


% % Logistic Loss
% n = length(Y);
% grad = zeros(size(W,2), 1);
% 
% for i=1:n
%     w = W(i,:);
%     y = Y(i);
%     
%     v = exp(-y* (w*x));
%     % grad = grad + (-y*v) /(1 + v) *(w');
%     grad = grad + (-y) /(1 + 1/v) *(w');
% end
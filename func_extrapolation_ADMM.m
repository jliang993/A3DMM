function e = func_extrapolation_ADMM(X, para)
% simple LP

type = para.type;

p = para.p; % numer of difference vectors
q = p - 2;

V = diff(X, 1,2);
V_k = V(:, 2:end);
v_k = V(:, end);
V_km1 = V(:, 1:end-1);

%%%% compute coefficient
VtV = (V_km1')*V_km1;
c = VtV \ ( (V_km1') * v_k );
% c = pinv(V_km1)* v_k;

%%%% iteration matrix
C = diag(ones(q-1,1), -1);
C(:, end) = c;

if isnan(sum(c))
    e = v_k;
else
    
    rho  = max(abs(eig(C)));
    % if rho>=1; disp(rho); end
    
    switch type
        case 'LP'
            s = para.s;
            if rho<1
                pC = C;
                for i=1:s+1; pC = pC*C; end
                S = (C - pC) / (eye(q) - C);
            else
                S = 0* C;
            end
            
        case 'LPinf'
            if rho<1
                S = C / (eye(q) - C);
                % S = eye(q) / (eye(q) - C) - C;
            else
                S = 0* C;
            end
            
    end
    
    e = V_k* S(:, end);
    
end


% EoF

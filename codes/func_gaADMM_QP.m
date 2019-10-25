function [lam,x,y, its, ek, dk, tk] = func_gaADMM_QP(para, xsol)

fprintf(sprintf('performing %s...\n', para.mname));
itsprint(sprintf('        step %08d: residual = %.3e...', 1,1), 1);

%%%%
if(~isfield(para,'p'))
    para.p = 4;
    para.gap = para.p + 1;
end

if(~isfield(para,'DoExtrapolation'))
    para.DoExtrapolation = 0;
end

if(~isfield(para,'type'))
    para.type = 'xxx';
end

if(~isfield(para,'SG'))
    safeguard_tol = 1e5;
else
    safeguard_tol = para.SG;
end

% dimension of problem, step-size
n = para.n;
gamma = para.gamma;
afun = para.afun;

q = para.p;
Q = para.Q;
Id = eye(n);
invM = Id /(gamma*Id + Q);

R = para.R;
L = para.L;

% stop cnd, max iteration
tol = para.tol;
maxits = para.maxits;

% past points
p = para.p; % # of past points
Z = zeros(prod(n), p);

%%%
dk = zeros(maxits, 1); % distance
ek = zeros(maxits, 1); % residual
tk = zeros(maxits, 1); % angle

% initial point and difference vector
x = zeros(prod(n), 1);
y = zeros(prod(n), 1);
lam = zeros(prod(n), 1);
z = lam + gamma*x;
zbar = z;

v = z(:) - z(:);

its = 1;
while(its<maxits)
    
    z_km1 = z;
    v_km1 = v;
    
    %%% y update
    y = max( min(zbar/gamma, R), L );
    
    %%% dual update
    lam = zbar -  gamma* y;
    
    %%% x update
    x = invM*(zbar - 2*lam - q);
    
    %%% z update
    z = lam + gamma*x;
    
    a = afun(its);
    zbar = z + a*(z - z_km1);
    
    
    v = z(:) - z_km1(:);
    
    Z = [Z(:, 2:end), z(:)];
    
    %%%%%% update states
    if 1
        res = norm(v(:));
        ek(its) = res;
        
        dk(its) = norm(x(:)-xsol(:));
                
        tk(its) = (v')*(v_km1) /norm(v)/norm(v_km1);
    end
    
    %%%%%% restarting, LP, MPE
    if (mod(its,para.gap)==0) && para.DoExtrapolation
        
        e_z = func_extrapolation_ADMM(Z, para);
        
        if para.SafeGuard
            coeff_z = min(1, safeguard_tol/(its^1.1*norm(e_z)));
        else
            coeff_z = 1;
        end
        
        z = z + coeff_z* e_z;
        zbar = z;
        
    end
    
    %%%%%% print information
    if mod(its,1e3)==0; itsprint(sprintf('      step %08d: residual = %.3e...', its,res), its); end
    
    %%%%%% stop?
    if (res<tol)||(res>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');


dk = dk(1:its-1);
ek = ek(1:its-1);
tk = tk(1:its-1);


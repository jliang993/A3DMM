function [psi,x,y, its, ek, dk_x, dk_z, tk] = func_gaADMM(para, xsol,zsol)

fprintf(sprintf('performing %s...\n', para.mname));
itsprint(sprintf('        step %08d: residual = %.3e...', 1,1), 1);

%%%%
if(~isfield(para,'p'))
    para.p = 2;
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
mu = para.mu;
gamma = para.gamma;
tau = mu /gamma;

proxR = para.proxR;
proxJ = para.proxJ;

% stop cnd, max iteration
tol = para.tol;
maxits = para.maxits;

afun = para.afun;

% past points
p = para.p; % # of past points
Z = zeros(prod(n), p);

%%%
dk_x = zeros(maxits, 1); % distance of x-x^star
dk_z = zeros(maxits, 1); % distance z-z^star
ek = zeros(maxits, 1); % residual
tk = zeros(maxits, 1); % angle

% initial point and difference vector
x = ones(prod(n), 1);
y = ones(prod(n), 1);
psi = ones(prod(n), 1);
z = psi + gamma*x;
zbar = z;

v = z(:) - z(:);

its = 1;
while(its<maxits)
    
    v_km1 = v;
    z_km1 = z;
    
    %%% y update
    y = proxJ(zbar/gamma, gamma);
    
    %%% dual update
    psi = zbar -  gamma* y;
    
    %%% x update
    x = proxR((zbar - 2*psi)/gamma, tau);
    
    %%% z update
    z = psi + gamma*x;
    
    %%% inertial step
    v = z - z_km1;
    
    a = afun(its);
    zbar = z + a* v;
    
    Z = [Z(:, 2:end), z(:)];
    
    %%%%%% update states
    res = norm(v(:));
    ek(its) = res;
    
    dk_x(its) = norm(x(:)-xsol(:));
    dk_z(its) = norm(z(:)-zsol(:));
    
    tk(its) = (v')*(v_km1) /norm(v)/norm(v_km1);
    
    %%%%%% Linear prediction
    if (mod(its,para.gap)==0) && para.DoExtrapolation
        
        e = func_extrapolation_ADMM(Z, para);
        
        if para.SafeGuard
            coeff = min(1, safeguard_tol/(its^1.1*norm(e)));
        else
            coeff = 1;
        end
        
        z = z + coeff* e;
        % y = proxR(z/gamma, tau);
        % psi = z -  gamma* y;
        zbar = z;
        
    end
    
    %%%%%% print information
    if mod(its,1e2)==0; itsprint(sprintf('      step %08d: residual = %.3e...', its,res), its); end
    
    %%%%%% stop?
    if (res<tol)||(res>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

dk_x = dk_x(1:its);
dk_z = dk_z(1:its);
ek = ek(1:its-1);
tk = tk(1:its-1);

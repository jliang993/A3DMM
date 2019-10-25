function [lamH,lamV,x,yH,yV, its, ek, dk, tk, qk] = func_gaADMM_INP(para, u, xsol)

verbose = para.verbose;

fprintf(sprintf('performing %s...\n', para.mname));
if verbose; itsprint(sprintf('        step %08d: residual = %.3e...', 1,1), 1); end

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
N = prod(n);

f = para.f;
S = para.S;
pS = para.pS;
ipS = 1 - pS;

gamma = para.gamma;
afun = para.afun;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = n(1); n2 = n(2);
dh = @(x) [ diff(x,1,2), zeros(n2,1) ];       % forward difference
dv = @(x) [ diff(x,1,1); zeros(1,n1) ];       % discrete y-derivative
dht = @(x) [ -x(:,1), -diff(x(:,1:n1-1),1,2), x(:,n1-1) ];      % transpose x-derivative
dvt = @(x) [ -x(1,:); -diff(x(1:n2-1,:),1,1); x(n2-1,:) ];      % transpose y-derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% stop cnd, max iteration
tol = para.tol;
maxits = para.maxits;

% past points
p = para.p; % # of past points
Z = zeros(2*prod(n), p);

%%%
dk = zeros(maxits, 1); % distance
ek = zeros(maxits, 1); % residual
tk = zeros(maxits, 1); % angle
qk = zeros(maxits, 1); % angle

% initial point and difference vector
x = ones(n);
yH = zeros(n);
yV = zeros(n);
lamH = zeros(n);
lamV = zeros(n);
zH = lamH + gamma*dh(x);
zV = lamV + gamma*dv(x);
zbarH = zH;
zbarV = zV;

gradF = @(x, zH,zV, lamH,lamV) dht(dh(x) - (zH-2*lamH)/gamma) + dvt(dv(x) - (zV-2*lamV)/gamma);
projS = @(x) f + ipS .* x;

z = [zH(:); zV(:)];
v = z(:) - z(:);

its = 1;
while(its<maxits)
    
    zH_km1 = zH;
    zV_km1 = zV;
    v_km1 = v;
    
    %%% y update
    yH = wthresh( zbarH /gamma, 's', 1/gamma);
    yV = wthresh( zbarV /gamma, 's', 1/gamma);
    
    %%% dual update
    lamH = zbarH -  gamma* yH;
    lamV = zbarV -  gamma* yV;
    
    %%% x update
    for i=1:S
        x_old = x;
        ix = x + (i-1)/(i+2)*(x - x_old);
        x = projS( ix - 1/8* gradF(ix, zbarH,zbarV, lamH,lamV) );
    end
    
    %%% z update
    zH = lamH + gamma*dh(x);
    zV = lamV + gamma*dv(x);
    
    a = afun(its);
    zbarH = zH + a*(zH - zH_km1);
    zbarV = zV + a*(zV - zV_km1);
    
    v = [zH(:); zV(:)] - [zH_km1(:); zV_km1(:)];
    
    z = [zH(:); zV(:)];
    
    Z = [Z(:, 2:end), z(:)];
    
    %%%%%% update states
    res = norm(v(:));
    ek(its) = res;
    
    qk(its) = psnr(uint8(x(:)), uint8(u(:)));
    
    dk(its) = norm(x(:)-xsol(:));
    
    tk(its) = (v')*(v_km1) /norm(v)/norm(v_km1);
    
    %%%%%% restarting, LP, MPE
    if (mod(its,para.gap)==0) && para.DoExtrapolation
        
        e_z = func_extrapolation_ADMM(Z, para);
        
        if para.SafeGuard
            coeff_z = min(1, safeguard_tol/(its^1.1*norm(e_z)));
        else
            coeff_z = 1;
        end
        
        e_zH = reshape(e_z(1:N), n);
        e_zV = reshape(e_z(N+1:end), n);
        
        zH = zH + coeff_z* e_zH;
        zV = zV + coeff_z* e_zV;
        zbarH = zH;
        zbarV = zV;
        
    end
    
    %%%%%% print information
    if mod(its,1e2)==0 && verbose; itsprint(sprintf('      step %08d: residual = %.3e...', its,res), its); end
    
    %%%%%% stop?
    if (res<tol)||(res>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');


dk = dk(1:its-1);
ek = ek(1:its-1);
tk = tk(1:its-1);
qk = qk(1:its-1);



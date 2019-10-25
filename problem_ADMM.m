function [para] = problem_ADMM(P)

caseR = P(1:end-3);
caseJ = P(end-1:end);

switch caseR
    %%%%%%%%%%%%%%%%%%%%%%
    case 'l1'
        k = 128;
        oversampling = 5;
        
        n = 2048; % length of signal
        m = oversampling*k; % # of measurements
        
        % original x
        x_true = 10* randn(n,1);
        x_true = sign(x_true) .* min( max(abs(x_true), 4), 16);
        mask = proj_mask(x_true,k/n, 'r');
        x_true = x_true .* mask;
        
        para.proxR = @(x, t) wthresh(x, 's', t);
        
        %%%%%%%%%%%%%%%%%%%%%%
    case 'l12'
        k = 128;
        oversampling = 5;
        
        n = 2048; % length of signal
        m = oversampling*k; % # of measurements
        
        % blk size and # of blks
        B = 4;
        kB= k/B;
        kappa = B;
        
        % blocksparse = @(B,IB) randn(B,length(IB));
        normB = @(x,B) sqrt(sum(reshape(x,[B length(x)/B]).^2,1));
        % e = @(x,B) x./reshape(repmat(normB(x,B),[B 1]),[L 1]);
        
        % original x
        IB = randperm(n/B);
        IB = IB(1:kB);
        x_true = zeros(B,n/B);
        for ib=1:length(IB)
            x_true(:,IB(ib)) = SparseVector(B, kappa, 'GAUSSIAN', true);
        end
        x_true = 4* x_true(:);
        
        para.B = B;
        
        para.proxR = @(x, t) max(1-t./reshape(repmat(normB(x,B),[B 1]),[n 1]),0).*x;
        
        %%%%%%%%%%%%%%%%%%%%%%
    case 'linfty'
        n = 512; % length of signal
        m = n - 2; % # of measurements
        
        % original x
        x_true = randn(n,1);
        
        max_ = max(abs(x_true)) + 1;
        
        s = 8;
        Rnd = randperm(n);
        
        x_true(Rnd(1:s)) = sign(x_true(Rnd(1:s))) .*max_;
        
        para.proxR = @(x, t) x - perform_l1ball_projection(x, t);
        
        %%%%%%%%%%%%%%%%%%%%%%
    case 'tv'
        k = 32;
        
        oversampling = 5;
        n = 512; % length of signal
        m = oversampling*k; % # of measurements
        
        % original x
        x_true = 5*SparseVector(n, k, 'GAUSSIAN', true);
        x_true = cumsum(x_true);
        x_true = x_true - mean(x_true);
        
        para.proxR = @(x, t) prox_tv1D(x, t);
        
        %%%%%%%%%%%%%%%%%%%%%%
    case 'lnuclear'
        n = [64, 64]; % size of the mtx
        r = 4; % rank of the mtx
        
        % m = min(3*r*(sum(n) - r), 1024);
        m = 3*r*(sum(n) - r);
        
        % original x
        x_true = randn(n(1), r) * randn(r, n(2)); % matrix to be tested
        x_true = x_true(:);
        
        para.proxR = @(x, t) vecsvt(x, t, n);
        
        %%%%%%%%%%%%%%%%%%%%%%
end

% forward operator
K = randn(m, prod(n)) /sqrt(m);

% observation
f = K* x_true(:);

if strcmp(caseJ, 'bp')
    
    pinvK = pinv(K);
    para.proxJ = @(x, gamma) x(:) + pinvK*(f-K*x(:));
    
else
    
    f = f + 1e-2*randn(size(f));
    
end

para.f = f;
para.n = n;
para.K = K;

para.mu = 1;

para.x_true = x_true;


% EOF
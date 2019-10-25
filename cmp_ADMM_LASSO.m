clear all
% close all
clc;

warning off;
addpath('toolbox') ;
set(groot,'defaultLineLineWidth', 1.0);
%%
strF = {'australian', 'mushrooms', 'covtype', 'ijcnn1',...
    'phishing', 'aloi', 'cod-rna', 'splice',...
    'breast-cancer', 'german-numer'};

i_file = 10;

%%% load and scale data
dataname = strF{i_file};
class_name = [dataname, '_label.mat'];
feature_name = [dataname, '_sample.mat'];

load(['data/', class_name]);
load(['data/', feature_name]);

K = full(h);
f = l;

[m, n] = size(K)

% rescale the data
fprintf(sprintf('rescale data...\n'));
itsprint(sprintf('      column %06d...', 1), 1);
for j=1:size(K,2)
    K(:,j) = rescale(K(:,j), -1, 1);
    if mod(j,1e2)==0; itsprint(sprintf('      column %06d...', j), j); end
end
fprintf(sprintf('\nDONE!\n\n'));
%%
para.f = f;
para.n = n; % dim of the problem
para.K = K;

para.mu = 1;

para.tol = 5e-11; % if change to 1e-11, LP-inf may have big oscillation
para.maxits = 1e5;
%%
para.gamma = 2;

Ktf = (K')*f /para.gamma;
Mj = eye(prod(para.n)) / (eye(prod(para.n)) + (K')*K/para.gamma);
para.proxJ = @(x, gamma) Mj*(Ktf + x);

para.proxR = @(x, t) wthresh(x, 's', t);
%% find solution
para.mname = 'normal ADMM';
para.DoExtrapolation = 0;

para.afun = @(k) 0;

[lamsol,xsol,ysol] = func_gaADMM_LASSO(para, 0,0);
zsol = lamsol + para.gamma* xsol;

[lam1,x1,y1, its1, ek1, dk_x1,dk_z1, tk1] = func_gaADMM_LASSO(para, xsol, zsol);

its1

fprintf('\n');
%% inertial ADMM
para.afun = @(k) 0.3;

[lamsol,xsol,ysol] = func_gaADMM_LASSO(para, 0,0);
zsol = lamsol + para.gamma* xsol;

[lam2,x2,y2, its2, ek2, dk_x2,dk_z2, tk2] = func_gaADMM_LASSO(para, xsol, zsol);

its2

fprintf('\n');
%% Linear prediction on ADMM, finite step, no line search
para.afun = @(k) 0;

para.DoExtrapolation = 1;

para.SG = 1e5;
para.SafeGuard = 1;

%%%% c
para.mname = 'ADMM + LP, finite';
para.type = 'LP'; % acceleration type

para.p = 6; % number of points, coeff c \in R^{p-2}
para.gap = para.p + 1;

para.s = 1e2; % LP steps

[lamsol,xsol,ysol] = func_gaADMM_LASSO(para, 0,0);
zsol = lamsol + para.gamma* xsol;

[lam3,x3,y3, its3, ek3, dk_x3,dk_z3, tk3] = func_gaADMM_LASSO(para, xsol, zsol);

its3

% semilogy(dk3)

fprintf('\n');
%% Linear prediction on ADMM, INfinite step
para.mname = 'ADMM + LP, INfinite';
para.type = 'LPinf'; % acceleration type

para.SG = 1e5;
para.SafeGuard = 1;

[lamsol,xsol,ysol] = func_gaADMM_LASSO(para, 0,0);
zsol = lamsol + para.gamma* xsol;

[lam4,x4,y4, its4, ek4, dk_x4,dk_z4, tk4] = func_gaADMM_LASSO(para, xsol, zsol);

its4

fprintf('\n');
%% plot angle      
axesFontSize = 10;
labelFontSize = 10;
legendFontSize = 11;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

figure(102), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);


p = semilogy(1-tk1, 'k', 'linewidth', 1.0);
axis([1, its1, min(1-tk1),1]);

grid on;
ax = gca;
ax.GridLineStyle = '--';


% ytick = [1e-12, 1e-8, 1e-4, 1e-0, 1e4];
% set(gca, 'yTick', ytick);

set(gca,'FontSize', 7);

ylb = ylabel({'$1-\cos(\theta_k)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.065, 0]);


filename = ['results', filesep, sprintf('admm-lasso-%s-tk.pdf', dataname)];
print(filename, '-dpdf');
%% plot dk
linewidth = 1;

axesFontSize = 10;
labelFontSize = 10;
legendFontSize = 11;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(101), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

xorz = 'z';
if strcmp(xorz, 'x')
    p = semilogy(dk_x1, 'color', [1/2,1/2,1/2], 'linewidth', 1);
    hold on;
    pi = semilogy(dk_x2, 'b', 'linewidth', 1);
    p3 = semilogy(dk_x3, 'r', 'linewidth', 1);
    p4 = semilogy(dk_x4, 'k--', 'linewidth', 1);
else
    p = semilogy(dk_z1, 'color', [1/2,1/2,1/2], 'linewidth', 1);
    hold on;
    pi = semilogy(dk_z2, 'b', 'linewidth', 1);
    p3 = semilogy(dk_z3, 'r', 'linewidth', 1);
    p4 = semilogy(dk_z4, 'k--', 'linewidth', 1);
end


grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, its1*1.01, 1e-9, max(dk_z1)]);
ytick = [1e-12, 1e-8, 1e-4, 1e-0, 1e4];
set(gca, 'yTick', ytick);


if strcmp(xorz, 'x')
    ylb = ylabel({'$\|x_{k}-x^\star\|$'}, 'FontSize', labelFontSize,...
        'FontAngle', 'normal', 'Interpreter', 'latex');
else
    ylb = ylabel({'$\|z_{k}-z^\star\|$'}, 'FontSize', labelFontSize,...
        'FontAngle', 'normal', 'Interpreter', 'latex');
end
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);

% lg = legend([p, p3, p4, p2_, pi],...
%     'ADMM', 'ADMM+LP, finite step', 'ADMM+LP, INfinite step', 'restarting FISTA', 'iADMM');
lg = legend([p, pi, p3, p4],...
    'ADMM', 'inertial ADMM, $a_{k}=0.3$', 'A$^3$DMM, $s=100$', 'A$^3$DMM, $s=+\infty$');
set(lg,'FontSize', legendFontSize);
% set(lg, 'Location', 'best');
set(lg, 'Interpreter', 'latex');
legend('boxoff');
%%%%%%%%%%%%%%%%

filename = ['results', filesep, sprintf('admm-lasso-%s.pdf', dataname)];
print(filename, '-dpdf');


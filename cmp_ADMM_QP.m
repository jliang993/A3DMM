clear all;
% close all;
clc

warning off;
addpath('toolbox') ;
%%
n = 512;

m = ceil( n /2 );
q = randn(m, n);
Q = (q')*q;

L = 4*randn(n, 1);
R = L + 4*rand(n, 1);
 
para.n = n;
para.q = randn(n, 1);
para.Q = Q;
para.L = L;
para.R = R;

para.verbose = 1;
para.tol = 5e-11; % stopping criterion
para.maxits = 5e4; % max # of iteration
%%
para.gamma = 10;
%% ADMM
para.mname = 'normal ADMM';

para.afun = @(k) 0;

[lamsol,xsol,ysol] = func_gaADMM_QP(para, 0);
% zsol = lamsol + para.gamma* xsol;

[lam1,x1,y1, its1, ek1, dk1, tk1] = func_gaADMM_QP(para, xsol);

its1

fprintf('\n');
%% inertial ADMM
para.mname = 'inertial ADMM';

para.afun = @(k) 0.3;

[lam2,x2,y2, its2, ek2, dk2, tk2] = func_gaADMM_QP(para, xsol);

its2

fprintf('\n');
%% Linear prediction on ADMM, finite step, no line search
para.afun = @(k) 0.0;

para.DoExtrapolation = 1;

para.SG = 1e5;
para.SafeGuard = 1;


%%%% c
para.mname = 'ADMM + LP, finite';
para.type = 'LP'; % acceleration type

para.p = 6; % number of points, coeff c \in R^{p-2}
para.gap = para.p + 1;

para.s = 1e2; % LP steps

[lamsol,xsol,ysol] = func_gaADMM_QP(para, 0);

[lam3,x3,y3, its3, ek3, dk3, tk3] = func_gaADMM_QP(para, xsol);

its3

fprintf('\n');
%% Linear prediction on ADMM, INfinite step
para.mname = 'ADMM + LP, INfinite';
para.type = 'LPinf'; % acceleration type

[lamsol,xsol,ysol] = func_gaADMM_QP(para, 0);

[lam4,x4,y4, its4, ek4, dk4, tk4] = func_gaADMM_QP(para, xsol);

its4

fprintf('\n');
%% distance error ||z_{k}-z^\star||
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

p1d = semilogy(dk1, 'color', [1/2,1/2,1/2], 'LineWidth',linewidth);
hold on,

p2d = semilogy(dk2, 'b', 'LineWidth',linewidth);

p3d = semilogy(dk3, 'r', 'LineWidth',linewidth);

p4d = semilogy(dk4, 'k--', 'LineWidth',linewidth);


uistack(p1d, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, its2, 1e-7, 2*max(dk1)]);
ytick = [1e-8, 1e-4, 1e-0, 1e4];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x^\star\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p1d, p2d, p3d, p4d], ...
    'ADMM', 'inertial ADMM, $a_k=0.3$',...
    'A$^3$DMM, $s = 100$', 'A$^3$DMM, $s = +\infty$');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
% set(lg, 'Location', 'best');
legend('boxoff');

filename = ['results', filesep, sprintf('admm-qp-dk.pdf')];
print(filename, '-dpdf');
%% angle
linewidth = 1;

axesFontSize = 10;
labelFontSize = 10;
legendFontSize = 14;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

figure(102), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

p1d = semilogy(1-tk1, 'color', [1/5,1/5,1/5], 'LineWidth',1);
hold on,


grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, its1, 1e-8, 1]);
ytick = [1e-12, 1e-8, 1e-4, 1e-0, 1e4];
set(gca, 'yTick', ytick);

ylb = ylabel({'$1-\cos(\theta_k)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


% lg = legend([p1d], '$1-\cos(\theta_k)$');
% set(lg,'FontSize', legendFontSize);
% set(lg, 'Interpreter', 'latex');
% % set(lg, 'Location', 'best');
% legend('boxoff');

filename = ['results', filesep, sprintf('admm-qp-thetak.pdf')];
print(filename, '-dpdf');
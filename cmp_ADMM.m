clear all;
% close all;
clc

warning off;
addpath('toolbox') ;
%% BP-type and LASSO-type
Problem = {'l1-bp'; 'l12-bp'; 'linfty-bp'; 'lnuclear-bp'; 'tv-bp';...    % 1-5
    'l1-lr'; 'l12-lr'; 'linfty-lr'; 'lnuclear-lr'; 'tv-lr';};  % 6-10
% l1: l1-norm
% l12: l12-norm
% linfty: linfty-norm
% lnuclear: nuclear norm
% tv: total variation
% bp: Basis Pursuit
% lr: Linear regression

prb = 1;
problem = Problem{prb};
para = problem_ADMM(problem);
para.problem = problem;

para.tol = 1e-10;
para.maxits = 2e4;
%%
para.gamma = 1;

if contains(problem, 'lr')
    para.mu = 2;
    
    f = para.f;
    K = para.K;
    
    Ktf = (K')*f /para.gamma;
    Mj = eye(prod(para.n)) / (eye(prod(para.n)) + (K')*K/para.gamma);
    para.proxJ = @(x, gamma) Mj*(Ktf + x);
end
%% ADMM
para.mname = 'normal ADMM';

para.afun = @(k) 0;

[lamsol,xsol,ysol] = func_gaADMM(para, 0, 0);
zsol = lamsol + para.gamma* xsol;

[lam1,x1,y1, its1, ek1, dk_x1,dk_z1, tk1] = func_gaADMM(para, xsol,zsol);

its1

fprintf('\n');
%% inertial ADMM
para.mname = 'inertial ADMM';

para.afun = @(k) 0.3; %(k-1)/(k+3);

[lamsol,xsol,ysol] = func_gaADMM(para, 0, 0);
zsol = lamsol + para.gamma* xsol;

[lam2,x2,y2, its2, ek2, dk_x2, dk_z2, tk2] = func_gaADMM(para, xsol,zsol);

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

para.s = 100; % LP steps

[lamsol,xsol,ysol] = func_gaADMM(para, 0, 0);
zsol = lamsol + para.gamma* xsol;

[lam3,x3,y3, its3, ek3, dk_x3,dk_z3, tk3] = func_gaADMM(para, xsol,zsol);

its3

fprintf('\n');
%% Linear prediction on ADMM, INfinite step
para.mname = 'ADMM + LP, INfinite';
para.type = 'LPinf'; % acceleration type

[lamsol,xsol,ysol] = func_gaADMM(para, 0, 0);
zsol = lamsol + para.gamma* xsol;

[lam4,x4,y4, its4, ek4, dk_x4,dk_z4, tk4] = func_gaADMM(para, xsol,zsol);

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

xorz = 'x';
if strcmp(xorz, 'x')
    p1d = semilogy(dk_x1, 'color', [1/2,1/2,1/2], 'LineWidth',linewidth);
    hold on,
    
    p2d = semilogy(dk_x2, 'b', 'LineWidth',linewidth);
    p3d = semilogy(dk_x3, 'r', 'LineWidth',linewidth);
    p4d = semilogy(dk_x4, 'k--', 'LineWidth',linewidth);
else
    p1d = semilogy(dk_z1, 'color', [1/2,1/2,1/2], 'LineWidth',linewidth);
    hold on,
    
    p2d = semilogy(dk_z2, 'b', 'LineWidth',linewidth);
    p3d = semilogy(dk_z3, 'r', 'LineWidth',linewidth);
    p4d = semilogy(dk_z4, 'k--', 'LineWidth',linewidth);
end

uistack(p1d, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, its1, 1e-9, 2*max(dk_z1)]);
ytick = [1e-8, 1e-4, 1e-0, 1e4];
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


lg = legend([p1d, p2d, p3d, p4d], ...
    'ADMM', 'inertial ADMM, $a_3=0.3$',...
    'A$^3$DMM, $s = 100$', 'A$^3$DMM, $s = +\infty$');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
% set(lg, 'Location', 'best');
legend('boxoff');

filename = ['results', filesep, sprintf('admm-bp-dk-%s.pdf', problem)];
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

p1d = plot(tk1, 'color', [1/5,1/5,1/5], 'LineWidth',1);
hold on,

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, its1, tk1(end)* 0.8, 1]);

ylb = ylabel({'$\cos(\theta_k)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


% lg = legend([p1d], '$\log(1-\cos(\theta_k))$');
lg = legend([p1d], '${\\cos}(\theta_k)$');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
set(lg, 'Location', 'SouthEast');
legend('boxoff');

filename = ['results', filesep, sprintf('admm-bp-thetak-%s.pdf', problem)];
print(filename, '-dpdf');
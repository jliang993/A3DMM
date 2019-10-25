clear all
% close all
clc;

warning off;
addpath('toolbox') ;
set(groot,'defaultLineLineWidth', 1.5);
%%
problem = 'l1-lr'; % LASSO problem
% problem = 'l1-bp'; % Basis Pursuit problem
[para] = problem_ADMM(problem);

para.tol = 1e-12; % tolerance
para.maxits = 3e3; % max number of iteration

para.problem = problem;

para.mu = 1; % penality parameter for l1-norm

f = para.f; % observation
K = para.K; % linear operator
normKtK = norm(K*K');
%%
% gamma, augmented Lagrangian coefficient
para.gamma = normKtK /10; sorf = 'f'; % sorf: success or failure
para.gamma = normKtK + 0.1; sorf = 's';

if contains(problem, 'bp'); sorf = 'f'; end

if contains(problem, 'lr')
    Ktf = (K')*f /para.gamma;
    Mj = eye(prod(para.n)) / (eye(prod(para.n)) + (K')*K/para.gamma);
    para.proxJ = @(x, gamma) Mj*(Ktf + x);
end
%% find solution
para.mname = 'normal ADMM';

para.afun = @(k) 0;

[lamsol,xsol,ysol] = func_gaADMM(para, 0,0);
zsol = lamsol + para.gamma*xsol;

[psi,x,y, its, ek, dk, ~, tk] = func_gaADMM(para, xsol, zsol);

its

fprintf('\n');

% figure(111); plot(tk); pause(1);
%% inertial ADMM
para.mname = 'inertial ADMM - 1';

para.afun = @(k) 0.3;

[psi1,x1,y1, its1, ek1, dk1, ~, ~] = func_gaADMM(para, xsol, zsol);

its1

fprintf('\n');
%% inertial ADMM
para.mname = 'inertial ADMM - 2';

para.afun = @(k) 0.7;

[psi2,x2,y2, its2, ek2, dk2, ~, ~] = func_gaADMM(para, xsol, zsol);

its2

fprintf('\n');
%% FISTA-ADMM
para.mname = 'FISTA-ADMM';

para.afun = @(k) (k-1)/(k+3);

[psi3,x3,y3, its3, ek3, dk3, ~, ~] = func_gaADMM(para, xsol, zsol);

its3

fprintf('\n');
%% plot dk
linewidth = 1;

axesFontSize = 6;
labelFontSize = 10;
legendFontSize = 10;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(101), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

p = semilogy(dk, 'color',[1/2,1/2,1/2], 'linewidth', 1.0);

hold on;

p1 = semilogy(dk1, 'm', 'linewidth', 1.0);

p2 = semilogy(dk2, 'r', 'linewidth', 1.0);

p3 = semilogy(dk3, 'b', 'linewidth', 1.0);

% p4 = semilogy(dk4, 'k', 'linewidth', 1.0);

uistack(p3, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, its, 1e-10, 2*max(dk)]);
ytick = [1e-12, 1e-8, 1e-4, 1e-0, 1e4];
set(gca, 'yTick', ytick);

set(gca,'FontSize', 7);

ylb = ylabel({'$\|z_k-z^\star\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.065, 0]);


lg = legend([p, p1, p2, p3],...
    'ADMM', 'inertial ADMM, $a_k=0.3$', 'inertial ADMM, $a_k=0.7$',...
    'inertial ADMM, $a_k=\frac{k-1}{k+3}$');
set(lg,'FontSize', legendFontSize);
set(lg, 'Location', 'SouthWest');
if para.gamma >= normKtK && contains(problem, 'lr')
    set(lg, 'Location', 'NorthEast');
else
    set(lg, 'Location', 'SouthWest');
end
set(lg, 'Interpreter', 'latex');
legend('boxoff');

if strcmp(sorf, 'f') || contains(problem, 'bp')
    filename = ['results', filesep, sprintf('failure-iADMM-%s.pdf', problem)];
else
    filename = ['results', filesep, sprintf('success-iADMM-%s.pdf', problem)];
end
print(filename, '-dpdf');
%%
figure(102), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

if strcmp(sorf, 'f') || contains(problem, 'bp')
    p = plot(tk, 'k', 'linewidth', 1.0);
    axis([1, its, 3/4,1]);
    
    ylb = ylabel({'$\cos(\theta_k)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
else
    p = semilogy(1-tk, 'k', 'linewidth', 1.0);
    axis([1, its, 1e-10,1]);
    
    ylb = ylabel({'$1 - \cos(\theta_k)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
end

grid on;
ax = gca;
ax.GridLineStyle = '--';

set(gca,'FontSize', 7);

set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.065, 0]);

if strcmp(sorf, 'f') || contains(problem, 'bp')
    filename = ['results', filesep, sprintf('failure-iADMM-%s-tk.pdf', problem)];
else
    filename = ['results', filesep, sprintf('success-iADMM-%s-tk.pdf', problem)];
end
print(filename, '-dpdf');

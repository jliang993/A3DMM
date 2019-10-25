clear all;
% close all;
clc

warning off;
addpath('data') ;
addpath('toolbox') ;
%%
% name = 'barbara.png';
name = 'cameraman.png';
% name = 'lena.png';

u = double(imread(name));

% u = u(1:2:end, 1:2:end); 

para.n = size(u);

para.verbose = 1;
para.tol = 1e-9; % stopping criterion
para.maxits = 5e3; % max # of iteration

para.gamma = 1;

pS = proj_mask(u, 0.5, 'p');
f = u .* pS;

para.f = f;
para.pS = pS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n1, n2] = size(u);
dh = @(x) [ diff(x,1,2), zeros(n2,1) ];       % forward difference
dv = @(x) [ diff(x,1,1); zeros(1,n1) ];       % discrete y-derivative
dht = @(x) [ -x(:,1), -diff(x(:,1:n1-1),1,2), x(:,n1-1) ];      % transpose x-derivative
dvt = @(x) [ -x(1,:); -diff(x(1:n2-1,:),1,1); x(n2-1,:) ];      % transpose y-derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADMM
para.S = 10;
para.DoExtrapolation = 0;

para.mname = 'normal ADMM';

para.afun = @(k) 0;

para.maxits = 1e4;

[~,~,xsol, ~,~, ~, ~, ~, tk1, qk1] = func_gaADMM_INP(para, u, 0);
[lamH1,lamV1,x1,yH1,yV1, its1, ek1, dk1, ~] = func_gaADMM_INP(para, u, xsol);

para.maxits = 9;

[~,~,x1] = func_gaADMM_INP(para, u, xsol);

its1

fprintf('\n');
%% inertial ADMM
para.DoExtrapolation = 0;

para.mname = 'inertial ADMM';

para.afun = @(k) 0.3;

para.maxits = 1e4;

[~,~,xsol, ~,~, ~, ~, ~, ~, qk2] = func_gaADMM_INP(para, u, 0);
[lamH2,lamV2,x2,yH2,yV2, its2, ek2, dk2, tk2] = func_gaADMM_INP(para, u, xsol);

para.maxits = 9;

[~,~,x2] = func_gaADMM_INP(para, u, xsol);

its2

fprintf('\n');
%% Linear prediction on ADMM, finite step, no line search
para.afun = @(k) 0.0;

para.DoExtrapolation = 1;

para.Theta = 1 + 1e-3;

para.SG = 1e5;
para.SafeGuard = 1;

para.IfScale = 0;
para.IfInertial = 0;

%%%% c
para.mname = 'ADMM + LP, finite';
para.type = 'LP'; % acceleration type

para.p = 6; % number of points, coeff c \in R^{p-2}
para.gap = para.p + 1;

para.s = 1e2; % LP steps

para.maxits = 1e4;

[~,~,xsol, ~,~, ~, ~, ~, ~, qk3] = func_gaADMM_INP(para, u, 0);

[lamH3,lamV3,x3,yH3,yV3, its3, ek3, dk3, tk3] = func_gaADMM_INP(para, u, xsol);

para.maxits = 9;

[~,~,x3] = func_gaADMM_INP(para, u, xsol);

its3

fprintf('\n');
%% Linear prediction on ADMM, INfinite step
para.mname = 'ADMM + LP, INfinite';
para.type = 'LPinf'; % acceleration type

% para.p = 6;
% para.gap = para.p + 0;

para.SG = 1e5;
para.SafeGuard = 1;

para.maxits = 1e4;


[~,~,xsol, ~,~, ~, ~, ~, ~, qk4] = func_gaADMM_INP(para, u, 0);
[lamH4,lamV4,x4,yH4,yV4, its4, ek4, dk4, tk4] = func_gaADMM_INP(para, u, xsol);

para.maxits = 9;

[~,~,x4] = func_gaADMM_INP(para, u, xsol);

its4

fprintf('\n');
%% distance error ||z_{k}-z^\star||
linewidth = 1;

axesFontSize = 8;
labelFontSize = 8;
legendFontSize = 11;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(111), clf;
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

axis([1, its1, 1e-9, 2*max(dk1)]);
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

filename = ['results', filesep, sprintf('admm-inp-dk-%s.pdf', name(1:end-4))];
print(filename, '-dpdf');
%% psnr
linewidth = 1;

axesFontSize = 10;
labelFontSize = 10;
legendFontSize = 11;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(111), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);


p1d = plot(qk1, 'color', [1/2,1/2,1/2], 'LineWidth',linewidth);
hold on,

p2d = plot(qk2, 'b', 'LineWidth',linewidth);

p3d = plot(qk3, 'r', 'LineWidth',linewidth);

p4d = plot(qk4, 'k--', 'LineWidth',linewidth);

uistack(p1d, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, 1e2, 22, 27.5]);
ytick = 22:27;
set(gca, 'yTick', ytick);

ylb = ylabel({'PSNR'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p1d, p2d, p3d, p4d], ...
    'ADMM', 'inertial ADMM, $a_k=0.3$',...
    'A$^3$DMM, $s = 100$', 'A$^3$DMM, $s = +\infty$');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
set(lg, 'Location', 'SouthEast');
legend('boxoff');

filename = ['results', filesep, sprintf('admm-inp-psnr-%s.pdf', name(1:end-4))];
print(filename, '-dpdf');
%%
figure(222), clf; 
plot(qk1, 'r')
hold on;
plot(qk2, 'k')
plot(qk3, 'm')
plot(qk4, 'b--')

grid on;
%%
figure(333), clf;

aa = 0;
subplot(221), imgsc(x1-aa*u); 
subplot(222), imgsc(x2-aa*u);
subplot(223), imgsc(x3-aa*u);
subplot(224), imgsc(x4-aa*u);

if 1
    resolution = 300; % output resolution
    output_size = 300 *[8, 8]; % output size
    
    figure(101), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(x1);
    
    filename = ['results', filesep, sprintf('step-8-ADMM.pdf')];
    print(filename, '-dpdf');
    
    
    figure(102), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(x2);
    
    filename = ['results', filesep, sprintf('step-8-iADMM.pdf')];
    print(filename, '-dpdf');
    
    figure(103), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(x3);
    
    filename = ['results', filesep, sprintf('step-8-sADMM.pdf')];
    print(filename, '-dpdf');
    
    figure(104), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(x4);
    
    filename = ['results', filesep, sprintf('step-8-infADMM.pdf')];
    print(filename, '-dpdf');
end
%% angle
linewidth = 1;

axesFontSize = 10;
labelFontSize = 10;
legendFontSize = 11;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(112), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

p1d = plot(tk1, 'color', [1/5,1/5,1/5], 'LineWidth',1);
hold on,


grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, its1, 3/4, 1]);
% ytick = [1e-8, 1e-4, 1e-0, 1e4];
% set(gca, 'yTick', ytick);

ylb = ylabel({'$\cos(\theta_k)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


% lg = legend([p1d], '$\cos(\theta_k)$');
% set(lg,'FontSize', legendFontSize);
% set(lg, 'Interpreter', 'latex');
% % set(lg, 'Location', 'best');
% legend('boxoff');

filename = ['results', filesep, sprintf('admm-inp-thetak-%s.pdf', name(1:end-4))];
print(filename, '-dpdf');
%%
if 1
    resolution = 300; % output resolution
    output_size = 300 *[8, 8]; % output size
    
    figure(101), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(u);
    
    filename = ['results', filesep, sprintf('%s-original-img.pdf', name(1:end-4))];
    print(filename, '-dpdf');
    % filename = ['results', filesep, sprintf('original-img.png')];
    % print(filename, '-dpng');
    
    
    figure(102), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(f);
    
    filename = ['results', filesep, sprintf('%s-corrupted-img.pdf', name(1:end-4))];
    print(filename, '-dpdf');
    % filename = ['results', filesep, sprintf('corrupted-img.png')];
    % print(filename, '-dpng');
    
    figure(103), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(x4);
    
    filename = ['results', filesep, sprintf('%s-inpainted-img.pdf', name(1:end-4))];
    print(filename, '-dpdf');
    % filename = ['results', filesep, sprintf('inpainted-img.png')];
    % print(filename, '-dpng');
end
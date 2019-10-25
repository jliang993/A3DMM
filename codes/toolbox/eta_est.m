function eta = eta_est(ek, p1, p2)

flag = 0;
if nargin<2
    flag = 1;
end

if flag
    figure(110),
    semilogy(ek,'r', 'LineWidth',1.25);
    grid on;
    title('\eta estimation');
    
    [p1, ~] = ginput(1);
    [p2, ~] = ginput(1);
    close(110);
    
    p1 = floor(p1);
    p2 = floor(p2);
end

a = ek(p1);
b = ek(p2);
eta = (b/a) ^ (1/(p2-p1));
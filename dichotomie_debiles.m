
a = rand(1);
b = 100*(randn(1)-0.5);
c = 100*(randn(1)-0.5);
f = @(x) a*x.^4 - b*x + c;



x0 = 0;
xp = 1;
fx0 = f(x0);
fxp = f(xp);

ca_monte = false;
while ~ca_monte
    if (fxp > fx0)
        ca_monte = true;
    else
        fx0 = fxp;
        xp = 10*xp;
        fxp = f(xp);
    end
end


x0 = 0;
xn = -1;
fx0 = f(x0);
fxn = f(xn);

ca_monte = false;
while ~ca_monte
    if (fxn > fx0)
        ca_monte = true;
    else
        fx0 = fxn;
        xn = 10*xn;
        fxn = f(xn);
    end
end

closest = [1, 2, 3;
    1, 2, 3;
    2, 3, 4; 
    3, 4, 5;
    3, 4, 5];

grid = linspace(xn, xp, 5);
funct = [fxn, 0, 0, 0, fxp];
funct(2) = f(grid(2));
funct(3) = f(grid(3));
funct(4) = f(grid(4));

compx = grid;
compf = funct;

k = 1;

while grid(5) - grid(1) > 1e-1
    [m, idx] = min(funct);
    grid = linspace(grid(closest(idx, 1)), grid(closest(idx, 3)), 5);
    temp = funct(closest(idx, :));
    funct(1) = temp(1);
    funct(3) = temp(2);
    funct(5) = temp(3);
    funct(2) = f(grid(2));
    funct(4) = f(grid(4));
    fprintf('iter=%07d x=%15.13e f=%15.13e\n', k, grid(idx), m);
    k = k+1;
    
    
    
    compx = cat(2, compx, grid(2), grid(4));
    compf = cat(2, compf, funct(2), funct(4));
    
end

xx = linspace(xn, xp, 100);
figure(1), clf;
plot(xx, f(xx), 'Color', 'blue', 'LineWidth', 5);
hold on;
plot([grid(idx),grid(idx)], [-1e10, 1e10], 'Color', 'red', 'LineWidth', 5);
axis([xn, xp, min(m, min(f(xn), f(xp))), max(f(xn), f(xp))]);


[~, idx] = sort(compx);
plot(compx(idx), compf(idx), 'Color', 'black', 'LineWidth', 5);
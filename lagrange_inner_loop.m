function [ us, Us, obj ] = lagrange_inner_loop( Hi, A, k, V, params )
%LAGRANGE_INNER_LOOP Summary of this function goes here
%   Detailed explanation goes here
STEP = 2;

mup = 1;
[~, ~, fmu0] = f(0, Hi, A, k, V, params);
[usp, Usp, fmup] = f(mup, Hi, A, k, V, params);

ca_monte = false;
while ~ca_monte
    if (fmup > fmu0)
        ca_monte = true;
    else
        fmu0    = fmup;
        mup     = STEP * mup;
        [usp, Usp, fmup] = f(mup, Hi, A, k, V, params);
    end
end


mun = -1;
[~, ~, fmu0] = f(0, Hi, A, k, V, params);
[usn, Usn, fmun] = f(mun, Hi, A, k, V, params);

ca_monte = false;
while ~ca_monte
    if (fmun > fmu0)
        ca_monte = true;
    else
        fmu0    = fmun;
        mun     = STEP * mun;
        [usn, Usn, fmun] = f(mun, Hi, A, k, V, params);
    end
end

% building the look-up table for 3-nearest neighbors
closest = [1, 2, 3;
    1, 2, 3;
    2, 3, 4;
    3, 4, 5;
    3, 4, 5];


% initializing the grid and evaluations
grid = linspace(mun, mup, 5);
funct = [fmun, 0, 0, 0, fmup];
resu = {usn, 0, 0, 0, usp};
resU = {Usn, 0, 0, 0, Usp};

[resu{2}, resU{2}, funct(2)] = f(grid(2), Hi, A, k, V, params);
[resu{3}, resU{3}, funct(3)] = f(grid(3), Hi, A, k, V, params);
[resu{4}, resU{4}, funct(4)] = f(grid(4), Hi, A, k, V, params);

compx = grid;
compf = funct;

[~, idx] = min(funct);

iter = 1;

while grid(5) - grid(1) > 1e-4
    
    grid = linspace(grid(closest(idx, 1)), grid(closest(idx, 3)), 5);
    temp = funct(closest(idx, :));
    
    funct(1) = temp(1);
    funct(3) = temp(2);
    funct(5) = temp(3);
    
    [resu{2}, resU{2}, funct(2)] = f(grid(2), Hi, A, k, V, params);
    [resu{4}, resU{4}, funct(4)] = f(grid(4), Hi, A, k, V, params);
    
    compx = cat(2, compx, grid(2), grid(4));
    compf = cat(2, compf, funct(2), funct(4));
    
    [m, idx] = min(funct);
    
    iter = iter+1;
end

us = resu{idx};
Us = resU{idx};
obj = funct(idx);


end

function [us, Us, obj] = f(x, Hi, A, k, V, params)

switch params.relaxation
    case 'spectral'
        [us, Us, obj] = spectral_inner_loop(Hi - x*ones(size(Hi)), A);
    case 'graph-cut'
        [us, Us, obj] = graph_cut(Hi- x*ones(size(Hi)), A);
end

obj = obj + x * (2 * k - V);

end
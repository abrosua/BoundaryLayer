%% MAIN PROGRAM
clc; clear; close all
% Viscous-Inviscid Interaction Solver
% 1. Inviscid solver using Potential flow method
% 2. Boundary layer approximation method with:
%       2.1. Laminar: Karman-Pohlhausen
%       2.2. Transition: Cebeci and Smith (1974)
%       2.3. Turbulent: Temporarily neglected

% Input values
airfoil_dir = 'naca1408.txt';
c = 1; % [m] chord length

mu = 1.81206e-5; % dynamic viscosity
rho = 1.225;    % [kg/m^3] density
nu = mu/rho;

%% Input geometry
[xb, yb, m, mp1, U_inf, alpha] = readData(airfoil_dir);
num_panel = m;
q = 0.5*rho*U_inf^2;
Re = U_inf*c/nu;

%% Viscous-Inviscid Iteration
threshold = 1e-05;
err = 1000; % initial error (to avoid NaN)
iter = 0;
% Initialization
old_dels = zeros(num_panel, 1);
g = zeros(num_panel+1, 1);

while err >= threshold
    %% Calculate inviscid term using Potential Flow method
    % TEMPORARY SUBROUTINE
    [x, y, gamma, vtan, cp] = panelMethod(xb,yb,m,mp1,alpha,g);
    
    % The results are U_in, cp, cl
    Ue = abs(U_inf*vtan');

    %% Calculate Boundary Layer for both viscous and turbulent regime
    % Obtaining stagnation point
    [~, stag] = min(abs(Ue));
    up = stag:num_panel;
    low = stag:-1:1;
    
    % Initialization
    deltas = zeros(num_panel, 1);
    thetas = zeros(num_panel, 1);
    cf = zeros(num_panel, 1);
    
    % Boundary layer solver
    [deltas(up), thetas(up), cf(up), trans_u, sp_u] = boundaryLayer...
        (Ue(up), x(up)', y(up)', c, nu);        % UPPER airfoil
    [deltas(low), thetas(low), cf(low), trans_l, sp_l] = boundaryLayer...
        (Ue(low), x(low)', y(low)', c, nu);     % LOWER airfoil
    % Transition and separation points
    trans = [stag-1+trans_u stag+1-trans_l];
    sp = [stag-1+sp_u stag+1-sp_l];
    
    %% Iteration
    % Calculate boundary condition for inviscid solver
    g(1:end-1) = 0.03.*deltas;
    g(end) = 0.03.*deltas(m);
    
    % Looping criteria
    err = sum(abs((deltas - old_dels)./old_dels)*100);
    old_dels(:) = deltas;
    
    % Error plotting
    iter = iter + 1;
    num_iter(iter) = iter;
    err_iter(iter) = err;

    fprintf('iteration %d ------ error = %.2d\n', iter, err);
end

%% Post-Processing
% Plot error iteration
figure; grid on
plot(num_iter, err_iter, 'ro-', 'linewidth', 2);
xlabel('Number of iteration'); ylabel('Total error [%]');
title('Error convergence');

% Trailing edge error handling
cp(1) = 1; cp(end) = 1;
vtan(1) = 0; vtan(end) = 0;
% Calculate cl, cd and plot some graphs.
[cl, cd] = postPro(m, alpha, stag, [nan nan], sp, vtan, x, y, cp, cf, deltas);

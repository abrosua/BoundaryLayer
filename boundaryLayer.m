% Subroutine for boundary layer calculation using,
% 1. Approximation method for Laminar regime
% 2. Karman momentum-integral method for Turbulent regime
% NOTE:
%       *transition checked using Cebeci and Smith (1974) method
%           which is an improvement of Michel's method

function [deltas, thetas, cf, trans, sp] = boundaryLayer...
    (Ue, xp, yp, c, nu)
    %% Initialization
    % Momentum thickness calculation
    eps = 1e-05; trans = nan;
    M = length(xp);
    thetas = zeros(M, 1);
    
    % Calculate dUe/dx
    dUe = zeros(M, 1);
    for i = 1:length(dUe)
        if i ~= length(Ue)
            dUe(i) = (Ue(i+1) - Ue(i))/((xp(i+1) - xp(i))/c);
        else
            dUe(i) = (Ue(i) - Ue(i-1))/((xp(i) - xp(i-1))/c);
        end
    end
    
    % Calculate ds (for integration)
    rp = [xp, yp];
    ds = zeros(M-1, 1);
    for i = 1:M-1
        ds(i) = abs(norm(rp(i+1,:) - rp(i,:)));
    end

    %% LAMINAR Regime
    
    % Stagnation point
    thetas(1) = sqrt(0.075*nu/(abs(dUe(1))));
    % Thetas for the other points
    for i = 2:M
        integral = 0;
        xt = 0;
        for j = 2:i
            integral = integral + (Ue(j)^5 + Ue(j-1)^5)*ds(j-1)/2;
            xt = xt + abs(xp(j) - xp(j-1));
        end
        thetas(i) = sqrt(0.45*nu*integral/(Ue(i)^6));
        
        % Checking TRANSITION POINT
        LHS = abs(Ue(i))*thetas(i)/nu;
        k = abs(Ue(i))*xt/nu;
        %RHS = 2.9*(k^0.4); % Michel (1952)
        RHS = 1.174*(1 + (22400/k))*(k^0.46); % Cebeci and Smith (1974)
        temp = abs(LHS - RHS);
        if isnan(trans*i) % Assign the index for transition point!
            if temp < eps
                trans = i;
            end
        end
    end

    % Calculate lambda (Pressure gradient parameter)
    lambda = (thetas.^2).*dUe/nu;

    % Calculate tau_wall (wall shear stress) and deltas (disp. thickness)
    L = zeros(M,1);
    H = zeros(M,1);
    for i = 1:M
        z = 0.25 - lambda(i);
        if lambda(i) < 0.1 && lambda(i) > 0
            L(i) = 0.22 + 1.57*lambda(i) - 1.8*lambda(i)^2;
            H(i) = 2.61 - 3.75*lambda(i) + 5.24*lambda(i)^2;
        elseif lambda(i) <= 0 && lambda(i) > -0.1
            L(i) = 0.22 + 1.402*lambda(i) + (0.018*lambda(i)/...
                (lambda(i)+0.107));
            H(i) = 2.088 + 0.0731/(lambda(i)+0.14);
        elseif lambda(i) >= 0.1 && lambda(i) <= 0.25
            L(i) = L(i+1);
            H(i)= 2.0 + 4.14*z - 83.5*z^2 + 854*z^3 - 3337*z^4 + 4576*z^5;
        else
            L(i) = L(i-1);
            H(i) = H(i-1);
        end
    end
    
    deltas = H.*thetas;
    cf = 2*nu*L./(Ue.*thetas);
    cf(thetas == 0) = 0;
    
    %% TURBULENT Regime
    
    % Initialization
    kapa = 0.41; B = 5.0;  % Coles & Hirst (1968)
    a = @(PI) (2 + 3.179*PI + 1.5*PI^2)/(kapa*(1 + PI));    % eqn. 6-119a
    Re_theta = @(PI, lambd, H) (1 + PI)*exp(kapa*lambd - kapa*B - 2*PI)/...
        (kapa*H);   % eqn. 6-119b
    
    %% Separation point evaluation
    for i = 1:M
        if lambda(i) <= -0.09
            sp = i;
            break
        else
            sp = M;
        end
    end
end
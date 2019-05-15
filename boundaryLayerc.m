% Subroutine for boundary layer calculation using,
% 1. Approximation method for Laminar regime
% 2. Karman momentum-integral method for Turbulent regime
% NOTE:
%       *transition checked using Cebeci and Smith (1974) method
%           which is an improvement of Michel's method

classdef boundaryLayerc
    methods(Static)
        function [deltas, thetas, cf, trans, sp] = main(Ue, xp, yp, c, nu)
            %% Initialization
            % Momentum thickness calculation
            eps = 1e-05; trans = nan;
            thetas = zeros(M, 1); deltas = zeros(M, 1); cf = zeros(M, 1);
            M = length(xp);
            
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

            %% Boundary layer iteration

            % Stagnation point
            thetas(1) = sqrt(0.075*nu/(abs(dUe(1))));
            lambda = (thetas(1).^2).*dUe(1)/nu;
            deltas(1) = H.*thetas(1);
            cf(1) = 2*nu*L./(Ue(1).*thetas(1));
            % Thetas for the other points
            for i = 2:M
                if isnan(trans*i)   % LAMINAR
                    integral = 0; xt = 0; prev = [0 0]; % Initialization
                    for j = 2:i
                        integral = integral + (Ue(j)^5 + Ue(j-1)^5)*...
                            ds(j-1)/2;
                        xt = xt + abs(xp(j) - xp(j-1));
                    end
                    thetas(i) = sqrt(0.45*nu*integral/(Ue(i)^6));

                    % Calculate lambda (Pressure gradient parameter)
                    lambda = (thetas(i).^2).*dUe(i)/nu;
                    [L, H] = laminar(lambda, prev);
                    prev = [L H];

                    % Calculate tau_wall (wall shear stress) and
                    %   deltas (disp. thickness)
                    deltas(i) = H.*thetas(i);
                    if thetas(i) ~= 0
                        cf(i) = 2*nu*L./(Ue(i).*thetas(i));
                    else
                        cf(i) = 0;
                    end
                    
                    % Checking TRANSITION POINT
                    test = transition(thetas(i), xt, Ue(i), nu);
                    if test < eps
                        trans = i;
                    end
                    
                else    % TURBULENT
                    
                end
            end

            %% Separation point evaluation
            for i = 1:M
                if lambda <= -0.09
                    sp = i;
                else
                    sp = M;
                end
            end            
        end
        function [L, H] = laminar(lambda, temp)
            %% Laminar subroutine
            z = 0.25 - lambda;
            if lambda < 0.1 && lambda > 0
                L = 0.22 + 1.57*lambda - 1.8*lambda^2;
                H = 2.61 - 3.75*lambda + 5.24*lambda^2;
            elseif lambda <= 0 && lambda > -0.1
                L = 0.22 + 1.402*lambda + (0.018*lambda/...
                    (lambda + 0.107));
                H = 2.088 + 0.0731/(lambda + 0.14);
            elseif lambda >= 0.1 && lambda <= 0.25
                L = temp(1);
                H= 2.0 + 4.14*z - 83.5*z^2 + 854*z^3 - 3337*z^4 + 4576*z^5;
            else
                L = temp(1);
                H = temp(2);
            end
        end
        function tes = turbulent(x)
            % Initialization
            kapa = 0.41; B = 5.0;  % Coles & Hirst (1968)
            a = @(PI) (2 + 3.179*PI + 1.5*PI^2)/...
                (kapa*(1 + PI));	% eqn. 6-119a
            Re_theta = @(PI, lambd, H) (1 + PI)*...
                exp(kapa*lambd - kapa*B - 2*PI)/(kapa*H);   % eqn. 6-119b
        end
        function test = transition(thetas, xt, Ue, nu)
            LHS = abs(Ue)*thetas/nu;
            k = abs(Ue)*xt/nu;
            %RHS = 2.9*(k^0.4); % Michel (1952)
            RHS = 1.174*(1 + (22400/k))*(k^0.46); % Cebeci and Smith (1974)
            test = abs(LHS - RHS);        
        end
    end
end
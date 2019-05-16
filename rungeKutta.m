% 4th Order Runge-Kutta to integrate dtheta/dx

function theta = rungeKutta(theta0, x0, xStop, cf, H, Ue, dUe)
    ft = @(t) cf/2 - (2+H)*t*dUe/Ue;    % Function handle for f(x, theta)

    N = 2;
    theta = theta0; thetap = theta0;
    
    eps = 1e-03; temp = 1;
    while temp > eps
        h = (xStop - x0)/N;
        for i = 1:N
            k1 = h*ft(theta);
            k2 = h*ft(theta + k1/2);
            k3 = h*ft(theta + k2/2);
            k4 = h*ft(theta + k3);
            theta = theta + (k1 + 2*k2 + 2*k3 + k4)/6;
        end
        temp = abs(theta - thetap);
        thetap = theta;
    end

end
% Subroutine for post-processing, such as:
%   1. Lift and drag coefficient calculation.
%   2. Pressure coefficient plotting
%   3. Velocity distribution plotting
%   4. Boundary layer plotting (displacement thickness)

function [cl, cd] = postPro(m, alpha, stag, trans, sp, vtan,...
    x, y, cp, cf, delta)
    %% Calculate lift and drag coefficient
    
    e = 1.7;
    cf = e*cf; delta = e*delta;
    cp_u = 0; cp_l = 0;
    cd_lam = 0; cd_turb = 0;
    % Upper airfoil
    for i = stag:m-1
        cp_u = cp_u + (cp(i+1) + cp(i))*(x(i+1)-x(i))/2;
        if i < trans(1)
            cd_lam = cd_lam + (cf(i+1) + cf(i))*(x(i+1)-x(i))/2;
        else
            cd_turb = cd_turb + (cf(i+1) + cf(i))*(x(i+1)-x(i))/2;
        end
    end
    % Lower airfoil
    for i = stag:-1:2
        cp_l = cp_l + (cp(i-1) + cp(i))*(x(i-1)-x(i))/2;
        if i > trans(2)
            cd_lam = cd_lam + (cf(i-1) + cf(i))*(x(i-1)-x(i))/2;
        else
            cd_turb = cd_turb + (cf(i-1) + cf(i))*(x(i-1)-x(i))/2;
        end
    end
    cl = (cp_l - cp_u)*cos(alpha*pi/180);
    cd = (cd_lam + cd_turb);
    % Displaying the results
    disp('Results: ');
    fprintf('Cl = %.4f\n', cl);
    fprintf('Cd = %.4f\n', cd);

    %% Plot pressure coefficient and velocity distribution

    figure; hold on; grid on
    plot(x(stag:m), -cp(stag:m), 'b-');
    plot(x(1:stag), -cp(1:stag), 'r-');
    axis([0 1 min(-cp) max(-cp)]);
    xlabel('x/c'); ylabel('-c_p');
    title('Coefficient of pressure distribution');
    legend('Upper', 'Lower');
    hold off

    figure; hold on; grid on
    plot(x(stag:m), abs(vtan(stag:m)), 'b-');
    plot(x(1:stag), abs(vtan(1:stag)), 'r-');
    axis([0 1 0 max(vtan)]);
    xlabel('x/c'); ylabel('v/U');
    title('Velocity distribution');
    legend('Upper', 'Lower');
    hold off

    % Plot shear at wall
    figure; hold on; grid on
    plot(x(stag:sp(1)), cf(stag:sp(1)), 'b-');
    plot(x(sp(2):stag), cf(sp(2):stag), 'r-');
    axis([0 1 0 0.05]);
    xlabel('x/c'); ylabel('cf');
    title('Coefficient of friction distribution');
    legend('Upper', 'Lower');
    hold off

    %% Plot airfoil with boundary layer
    %boundary layer thickness
    yBL = zeros(m, 1);
    yBL(stag+1:end) = y(stag+1:end)' + delta(stag+1:end);
    yBL(1:stag-1) = y(1:stag-1)' - delta(1:stag-1);
    yBL(stag) = y(stag) + delta(stag)*sin(alpha+pi);
    xBL(stag) = x(stag) - delta(stag)*cos(alpha);

    figure; grid on; hold on; axis equal
    plot(x, y, 'k');
    %plot(x(trans(1):trans(2)), yBL(trans(1):trans(2)), 'r');
    %plot(x(trans(1)), yBL(trans(1):trans(2)), 'r');
    plot(x(sp(2)+1:sp(1)-1), yBL(sp(2)+1:sp(1)-1), 'r');
    plot(xBL(stag), yBL(stag), 'b*');
    plot(x(trans(1)), yBL(trans(1)), 'bo');
    plot(x(trans(2)), yBL(trans(2)), 'bo');
    legend('Airfoil', 'Boundary Layer', 'Stagnation Point',...
        'Transition Points');
    xlabel('x/c'); ylabel('y');
    title('Boundary Layer at Airfoil');

end

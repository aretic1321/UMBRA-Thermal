function f = calc_true_anomaly(t,mu,a,e)

    % Define convergence tolerance
    tol = 1e-6*pi/180;

    % Calculate orbit frequency
    n = sqrt(mu/a^3);

    % Calculate mean anomaly
    M = n*t;

    % Numerically estimate eccentric anomaly
    E1 = -inf;
    E2 = M;
    while abs(E2-E1) >= tol
        E1 = E2;
        E2 = E1 - (E1-e*sin(E1)-M)/(1-e*cos(E1));
    end

    % Calculate true anomaly
    f = 2*atan2(sqrt(1+e)*sin(E2/2), sqrt(1-e)*cos(E2/2));

end
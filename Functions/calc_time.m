function t = calc_time(f,mu,a,e)

    % Calculate eccentric anomaly
    E = 2*atan2(sqrt(1-e)*sin(f/2), sqrt(1+e)*cos(f/2));

    % Calculate mean anomaly
    M = E - e*sin(E);

    % Calculate orbit frequency
    n = sqrt(mu/a^3);

    % Calculate time after periapse
    t = M/n;

end
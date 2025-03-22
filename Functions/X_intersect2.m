function [x_sol, y_sol, z_sol] =...
    X_intersect2(d, R_s, gamma, phi, omega, psi)
    % X_intersect2 finds if and where there is an intersection


    % recovers theta, R, and h from the given values
    theta = asind(R_s/d);
    R = R_s*cosd(theta);
    h = R*cotd(theta);

    tol = 1e-14;
    % determines where an intersection occurs, if at all
    if omega == 0 || omega == 180 ||...
            (gamma==0 && phi==omega)...
            || ((gamma==180 || gamma==-180) && phi==180-omega)...
            || (phi - psi >= 90 - theta)...
            || (omega < 90 && 90 - omega >= theta)...
            || (omega > 90 && 90 - (180 - omega) >= theta)

        % guaranteed situations to have no solution
        x_sol = [];
        y_sol = [];
        z_sol = [];
    else
        if gamma == 0
            %{
            if phi == 90
                if psi == 90
                    x_sol = 0;
                else
                    x_sol = -R_s*cosd(psi);
                end

                z_sol = -x_sol*tand(omega);
            else
                if omega == 90
                    x_sol = 0;
                else
                    x_sol = (-d*cosd(omega)+R_s*cosd(omega)*cosd(psi)/cosd(phi))/...
                        (sind(omega)-tand(phi)*cosd(omega));
                end
                z_sol = d-R_s*cosd(psi)/cosd(phi)-x_sol*tand(phi);
            end
            %}
            if phi == 90
                x_sol = -R_s*cosd(psi);

                z_sol = -x_sol*tand(omega);
            else
                x_sol = (-d*cosd(omega)+R_s*cosd(omega)*cosd(psi)/cosd(phi))/...
                    (sind(omega)-tand(phi)*cosd(omega));

                z_sol = d-R_s*cosd(psi)/cosd(phi)-x_sol*tand(phi);
            end

            y_sol = zeros([2, 1]);
            y_sol(1) = sqrt(R_s^2-x_sol^2-(z_sol-d)^2);
            y_sol(2) = -y_sol(1);

            x_sol = repmat(x_sol, 2, 1);
            z_sol = repmat(z_sol, 2, 1);
        elseif gamma == 180 || gamma == -180
            %{
            if phi == 90
                if psi == 90
                    x_sol = 0;
                else
                    x_sol = -R_s*cosd(psi);
                end
                    
                z_sol = x_sol*tand(omega);
            else
                if omega == 90
                    x_sol = 0;
                else
                    x_sol = (-d*cosd(omega)+R_s*cosd(omega)*cosd(psi)/cosd(phi))/...
                        (-sind(omega)-tand(phi)*cosd(omega));
                end

                z_sol = d-R_s*cosd(psi)/cosd(phi)-x_sol*tand(phi);
            end
            %}
            if phi == 90
                x_sol = -R_s*cosd(psi);
                    
                z_sol = x_sol*tand(omega);
            else
                x_sol = (-d*cosd(omega)+R_s*cosd(omega)*cosd(psi)/cosd(phi))/...
                        (-sind(omega)-tand(phi)*cosd(omega));

                z_sol = d-R_s*cosd(psi)/cosd(phi)-x_sol*tand(phi);
            end
            y_sol = zeros([2, 1]);
            y_sol(1) = sqrt(R_s^2-x_sol^2-(z_sol-d)^2);
            y_sol(2) = -y_sol(1);

            x_sol = repmat(x_sol, 2, 1);
            z_sol = repmat(z_sol, 2, 1);
        else
            if phi == 90
               x_sol = -R_s*cosd(psi);
               %{
               if omega == 90
                   a = 1;
               else
                   a = cotd(omega)^2/sind(gamma)^2+1;
               end
               if omega == 90 || psi == 90
                   b = -2*d;
               else
                   b = -2*R_s*cosd(psi)*cotd(gamma)*cotd(omega)/sind(gamma)-2*d;
               end
               if psi == 90
                   c = d^2-R_s^2;
               else
                   c = R_s^2*cosd(psi)^2+R_s^2*cosd(psi)^2*cotd(gamma)^2+...
                       d^2-R_s^2;
               end
               %}
               
               a = cotd(omega)^2/sind(gamma)^2+1;
               b = -2*R_s*cosd(psi)*cotd(gamma)*cotd(omega)/sind(gamma)-2*d;
               c = R_s^2*cosd(psi)^2+R_s^2*cosd(psi)^2*cotd(gamma)^2+...
                       d^2-R_s^2;

               z_sol = [(-b+sqrt(b^2-4*a*c))/(2*a);...
                   (-b-sqrt(b^2-4*a*c))/(2*a)];

               y_sol = R_s*cosd(psi)*cotd(gamma)-z_sol*cotd(omega)/sind(gamma);

               x_sol = repmat(x_sol, 2, 1);
            else
                %{
                if phi == 0 || phi == 180
                    a = 1+cotd(gamma)^2;
                elseif omega == 90
                    a = 1+tand(phi)^2+cotd(gamma)^2;
                else
                    a=1+tand(phi)^2+(cotd(gamma)-cotd(omega)*...
                        tand(phi)/sind(gamma))^2;
                end

                b=(2*R_s*cosd(psi)*tand(phi))/cosd(phi)+(2*cosd(omega)*...
                    (d-(R_s*cosd(psi))/cosd(phi))*(cosd(gamma)*sind(omega)-...
                    cosd(omega)*tand(phi)))/(sind(gamma)^2*sind(omega)^2);

                if omega == 90 && psi == phi
                    c = 0;
                else
                    c=(R_s^2*cosd(psi)^2)/cosd(phi)^2-R_s^2+(cosd(omega)^2*....
                        (d-(R_s*cosd(psi))/cosd(phi))^2)/...
                        (sind(gamma)^2*sind(omega)^2);
                end
                %}
                a=1+tand(phi)^2+(cotd(gamma)-cotd(omega)*...
                        tand(phi)/sind(gamma))^2;
                b=(2*R_s*cosd(psi)*tand(phi))/cosd(phi)+(2*cosd(omega)*...
                    (d-(R_s*cosd(psi))/cosd(phi))*(cosd(gamma)*sind(omega)-...
                    cosd(omega)*tand(phi)))/(sind(gamma)^2*sind(omega)^2);
                c=(R_s^2*cosd(psi)^2)/cosd(phi)^2-R_s^2+(cosd(omega)^2*....
                    (d-(R_s*cosd(psi))/cosd(phi))^2)/...
                    (sind(gamma)^2*sind(omega)^2);

                x_sol = [(-b+sqrt(b^2-4*a*c))/(2*a);...
                   (-b-sqrt(b^2-4*a*c))/(2*a)];

                z_sol = d-R_s*cosd(psi)/cosd(phi)-x_sol*tand(phi);
                
                y_sol = (-x_sol*...
                    (sind(omega)*cosd(gamma)-tand(phi)*cosd(omega))-...
                    d*cosd(omega)+R_s*cosd(psi)*cosd(omega)/cosd(phi))/...
                    (sind(omega)*sind(gamma));
            end
        end

        % reduces the solution to fit conditions
        x_sol_temp = x_sol;
        y_sol_temp = y_sol;
        if any(isnan([x_sol, y_sol, z_sol]), 'all')
            x_sol = x_sol(~isnan(x_sol_temp) &...
                ~isnan(y_sol_temp) & ~isnan(z_sol));
            y_sol = y_sol(~isnan(x_sol_temp) &...
                ~isnan(y_sol_temp) & ~isnan(z_sol));
            z_sol = z_sol(~isnan(x_sol_temp) &...
                ~isnan(y_sol_temp) & ~isnan(z_sol));
        end
        x_sol_temp = x_sol;
        y_sol_temp = y_sol;
        if any(imag([x_sol, y_sol, z_sol]) ~= 0, 'all')
            x_sol = real(x_sol(abs(imag(x_sol_temp)) < tol &...
                abs(imag(y_sol_temp)) < tol &...
                abs(imag(z_sol)) < tol));
            y_sol = real(y_sol(abs(imag(x_sol_temp)) < tol &...
                abs(imag(y_sol_temp)) < tol &...
                abs(imag(z_sol)) < tol));
            z_sol = real(z_sol(abs(imag(x_sol_temp)) < tol &...
                abs(imag(y_sol_temp)) < tol &...
                abs(imag(z_sol)) < tol));
            if length(x_sol) == 2
                if abs(x_sol(2) - x_sol(1)) < tol &&...
                       abs(y_sol(2) - y_sol(1)) < tol &&...
                       abs(z_sol(2) - z_sol(1)) < tol
                    x_sol = x_sol(1);
                    y_sol = y_sol(1);
                    z_sol = z_sol(1);
                end
            end
        end
        x_sol = x_sol(z_sol >= d-R_s & z_sol < h);
        y_sol = y_sol(z_sol >= d-R_s & z_sol < h);
        z_sol = z_sol(z_sol >= d-R_s & z_sol < h);
        if length(x_sol) == 2
            % if the solutions are roughly equal, assume no
            % intersections sindce the integration across the
            % intersections is 0
            if abs(x_sol(2) - x_sol(1)) < tol &&...
                   abs(y_sol(2) - y_sol(1)) < tol &&...
                   abs(z_sol(2) - z_sol(1)) < tol
                x_sol = [];
                y_sol = [];
                z_sol = [];
            end
        end
    end
end

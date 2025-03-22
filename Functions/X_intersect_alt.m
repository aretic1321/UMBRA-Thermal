function [x_sol, y_sol, z_sol] =...
    X_intersect_alt(d, R_s, gamma, phi, omega, psi)
    % X_intersect finds if and where there is an intersection and creates
    % table enteries of generic solutions with the conditions for the
    % intersection to exist
    
    % recovers theta, R, and h from the given values
    theta = asin(R_s./d);
    R = R_s.*cos(theta);
    h = R.*cot(theta);

    % determines where an intersection occurs, if at all
    if omega == 0 || omega == pi ||...
            (gamma==0 && phi==omega)...
            || ((gamma==pi || gamma==-pi) && phi==pi-omega)
        % guaranteed situations to have no solution or a situation where
        % the intersection of the solution wouldn't yield anything
        % different from a situation with no solution for near parameters
        x_sol = [];
        y_sol = [];
        z_sol = [];
    else
        
        g = sym(gamma);
        o = sym(omega);
        ps = sym(psi);
        ph = sym(phi);
        Rss = sym(R_s);
        ds = sym(d);
        t = asin(Rss./ds);
        Rsm = Rss.*cos(t);
        hs = Rsm.*cot(t);
        
        x = [];
        y = [];
        z = [];
        syms x y z
        
        eq1 = x.^2+y.^2+(z-ds).^2==Rss.^2;
        eq2 = x.*sin(ph)+z.*cos(ph)==...
            ds.*cos(ph)-Rss.*cos(ps);
        eq3 = x.*sin(o).*cos(g)+...
        y.*sin(o).*sin(g)+...
        z.*cos(o)==0;
        
        assumeAlso(z <= hs)
        [x_sol, y_sol, z_sol, pams, constrs] = solve([eq1, eq2, eq3],...
            [x, y, z], real=true, ReturnConditions=true, IgnoreProperties=false);
        
        if ~isempty(pams)
            assume(constrs)
            restrs = [z_sol>=ds-Rss, z_sol <= hs];
            [pa1, pa2, pa3] = solve(restrs, pams);
            x_sol = subs(x_sol, pams(1), pa1);
            y_sol = subs(y_sol, pams(2), pa2);
            z_sol = subs(z_sol, pams(3), pa3);
            assume(constrs, 'clear')
        end
        x_sol = double(x_sol);
        y_sol = double(y_sol);
        z_sol = double(z_sol);
        if (gamma == 0 || gamma == pi || gamma == -pi)
            % in the case where there is one solution for gamma==0, pi, or
            % -pi, the solution is switched to no solution because it would
            % produce the same result (theoretically), but it is easier to
            % work with no solution.
                if isscalar(y_sol)
                    x_sol = [];
                    y_sol = [];
                    z_sol = [];
                elseif length(y_sol)==2
                    if abs(y_sol(2)-y_sol(1)) < 1e-12
                        x_sol = [];
                        y_sol = [];
                        z_sol = [];
                    end
                end
        end
    end
end

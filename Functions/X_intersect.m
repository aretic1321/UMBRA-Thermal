function [x_sol, y_sol, z_sol, X_sol_set] =...
    X_intersect(d, R_s, gamma, phi, omega, psi, X_sol_set)
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
        % guaranteed situations to have no solution
        x_sol = [];
        y_sol = [];
        z_sol = [];
    else
        clmntocmp = "";
        % determines gamma, sets up the column name with gamma
        if gamma == 0
            clmntocmp = clmntocmp + "g=0";
    
            % special case can occur at phi=omega
            if phi == omega
                clmntocmp = clmntocmp + "ph=o";
                %{
                if d.*cos(phi) == R_s.*cos(psi)
                    % s=special case
                    clmntocmp = clmntocmp + "s";
                end % other case is no solution
                %}
            else
                clmntocmp = clmntocmp + "ph~=o";
            end
        elseif gamma == pi || gamma == -pi
            clmntocmp = clmntocmp + "g=pig=-pi";
    
            % special case can occur at phi=pi - omega
            if phi == pi - omega
                clmntocmp = clmntocmp + "ph=pi-o";
                %{
                if abs(d.*cos(phi)-R_s.*cos(psi))<1e-14
                    % s=special case
                    clmntocmp = clmntocmp + "s";
                end % other case is no solution
                %}
            else
                clmntocmp = clmntocmp + "ph~=pi-o";
            end
        else
            clmntocmp = clmntocmp + "g~=0g~=pig~=-pi";
        end
        % determines phi, sets up the column name with phi
        if phi == pi/2
            clmntocmp = clmntocmp + "ph=pi/2";
        elseif phi < pi/2
            clmntocmp = clmntocmp + "ph<pi/2";
        else %if phi > pi/2
            clmntocmp = clmntocmp + "ph>pi/2";
        end
        %{
        if phi == 0
            clmntocmp = clmntocmp + "ph=0";
        elseif phi == pi
            clmntocmp = clmntocmp + "ph=pi";
        elseif phi == pi/2
            clmntocmp = clmntocmp + "ph=pi/2";
        elseif phi < pi/2
            clmntocmp = clmntocmp + "ph<pi/2";
        else %if phi > pi/2
            clmntocmp = clmntocmp + "ph>pi/2";
        end
        %}
        % determines omega, sets up the column name with omega
        if omega == pi/2
            clmntocmp = clmntocmp + "o=pi/2";
        elseif omega < pi/2
            clmntocmp = clmntocmp + "o<pi/2";
        else %if omega > pi/2
            clmntocmp = clmntocmp + "o>pi/2";
        end
        % determines psi, sets up the column name with psi (if beyond pi/2)
        if psi > pi/2
            clmntocmp = clmntocmp + "ps>pi/2";
        end

        g = [];
        o = [];
        ps = [];
        ph = [];
        Rss = [];
        ds = [];
        t = [];
        hs = [];
        syms g o ps ph Rss ds t hs
        if any(strcmp(clmntocmp,...
            X_sol_set.Properties.VariableNames))
            x_sol = X_sol_set{'x_sol', clmntocmp}{1};
            y_sol = X_sol_set{'y_sol', clmntocmp}{1};
            z_sol = X_sol_set{'z_sol', clmntocmp}{1};
            conx = X_sol_set{'conx', clmntocmp}{1};
            cony = X_sol_set{'cony', clmntocmp}{1};
            conz = X_sol_set{'conz', clmntocmp}{1};
        else
            x = [];
            y = [];
            z = [];
            syms x y z
            assume([x y z g o ps ph Rss ds t hs], 'real')
            assumeAlso(Rss>0)
            assumeAlso(ds>Rss)
            assumeAlso(z >= ds-Rss)
            assumeAlso(t==asin(Rss/ds))
            assumeAlso(hs == Rss.*cos(t).*cot(t))
            assumeAlso(z <= hs)
            if psi > pi/2
                assumeAlso(ps>pi/2 & ps <= pi)
            else
                assumeAlso(ps>0 & ps<=pi/2)
            end
            
            eq1 = x.^2+y.^2+(z-ds).^2==Rss.^2;
            eq2 = x.*sin(ph)+z.*cos(ph)==...
                ds.*cos(ph)-Rss.*cos(ps);
            eq3 = x.*sin(o).*cos(g)+...
            y.*sin(o).*sin(g)+...
            z.*cos(o)==0;
    
            assumeAlso(o>0 & o<pi)
            % determines gamma, sets up the column name with gamma
            if gamma == 0
                assumeAlso(g==0)    
                % special case can occur at phi=omega
                if phi == omega
                    assumeAlso(ph==o) % specific assumption
                    if d.*cos(phi) == R_s.*cos(psi)
                        assumeAlso(ds.*cos(ph) == Rss.*cos(ps))
                    end % other case is no solution
                else
                    % assumption to avoid special case
                    assumeAlso(ph~=pi - o)
                end
            elseif gamma == pi || gamma == -pi
                assumeAlso(g==pi | g==-pi)
        
                % special case can occur at phi=pi - omega
                if phi == pi - omega
                    assumeAlso(ph==pi - o) % specific assumption
                    %{
                    if abs(d.*cos(phi)-R_s.*cos(psi))<1e-14
                        assumeAlso(ds.*cos(ph) == Rss.*cos(ps))
                    end % other case is no solution
                    %}
                else
                    % assumption to avoid special case
                    assumeAlso(ph~=pi - o)
                end
            else
                assumeAlso(g~=0 & g~=pi & g~=-pi)
            end
            % determines phi, sets up the column name with phi
            if phi == pi/2
                assumeAlso(ph == pi/2)
            elseif phi < pi/2
                assumeAlso(ph>=0 & ph<pi/2)
            else %if phi > pi/2
                assumeAlso(ph>pi/2 & ph<=pi)
            end
            %{
            if phi == 0
                assumeAlso(ph == 0)
            elseif phi == pi
                assumeAlso(ph == pi)
            elseif phi == pi/2
                assumeAlso(ph == pi/2)
            elseif phi < pi/2
                assumeAlso(ph>0 & ph<pi/2)
            else %if phi > pi/2
                assumeAlso(ph>pi/2 & ph<pi)
            end
            %}
            % determines omega, sets up the column name with omega
            if omega == pi/2
                assumeAlso(o==pi/2)
            elseif omega < pi/2
                assumeAlso(o<pi/2)
            else %if omega > pi/2
                assumeAlso(o>pi/2)
            end
            
            eq2 = simplify(eq2);
            eq3 = simplify(eq3);
        
            % gamma==0, gamma==pi and gamma==-pi have special case 
            % scenarios where the plates made from eq2 and eq3 can be
            % parallel
            if gamma==0 || (gamma==pi || gamma==-pi)
                % special case can occur at gamma=0, phi=omega and
                % d.*cos(phi)=R_s.*cos(psi)
                % special case can occur at gamma=pi or -pi, phi=pi-omega 
                % and d.*cos(phi)=R_s.*cos(psi)
                %{
                if (gamma==0 && phi==omega) || ((gamma==pi ||...
                        gamma==-pi) && phi==pi-omega)
                    
                    if abs(d.*cos(phi)-R_s.*cos(psi))<1e-14
                        x_sol = solve(subs(eq3, z, hs), x);
                        [y_sol, ~, cony] = solve(...
                            subs(eq1, {z, x}, {hs, x_sol}), y,...
                            ReturnConditions=true);
                        conx = symtrue;
                        conz = symtrue;
                        
                        % repeat the x and z for the different y
                        % solutions
                        x_sol = repmat(x_sol, 2, 1);
                        z_sol = [hs; hs];
                        conx = repmat(conx, 2, 1);
                        conz = repmat(conz, 2, 1);
                        
                        cony=simplify(cony);
                        X_sol_set.(clmntocmp)=...
                            {x_sol; y_sol; z_sol; conx; cony; conz};
                    end
                else
                %}
                
                    if phi == pi/2
                        % solve for x
                        [x_sol, ~, conx] = solve(eq2, x,...
                                ReturnConditions=true);
                        % solve for z
                        [z_sol, ~, conz] = solve(...
                            subs(eq3, x, x_sol), z,...
                                ReturnConditions=true);
                        % solve for y
                        [y_sol, ~, cony] = solve(...
                            subs(eq1, {x, z}, {x_sol, z_sol}),...
                            y, ReturnConditions=true);
                        
                        if isempty(y_sol)
                            x_sol = [];
                            z_sol = [];
                        % repeat the x and y for the possibly different z
                        % solutions
                        elseif length(y_sol) == 2
                            if length(x_sol) < 2
                                x_sol = repmat(x_sol, 2, 1);
                            end
                            if length(conx) < 2
                                conx = repmat(conx, 2, 1);
                            end
                            if length(z_sol) < 2
                                z_sol = repmat(z_sol, 2, 1);
                            end
                            if length(conz) < 2
                                conz = repmat(conz, 2, 1);
                            end
                        end
                    else
                        % solve for z in terms of x
                        [z_sol, ~, conz] = solve(eq2, z,...
                                ReturnConditions=true);
                        % solve for x (maybe) in terms of y
                        [x_sol, ~, conx] = solve(...
                            subs(eq3, z, z_sol), x,...
                                ReturnConditions=true);
                        % gets z in terms of y (if possible)
                        z_sol = subs(z_sol, x, x_sol);
                        conz = subs(conz, x, x_sol);
                        % solve for y
                        [y_sol, ~, cony] = solve(...
                            subs(eq1, {z, x}, {z_sol, x_sol}),...
                            y, ReturnConditions=true);
                        z_sol = subs(z_sol, y, y_sol);
                        x_sol = subs(x_sol, y, y_sol);
                        %conz = simplify(subs(conz, y, y_sol));
                        %conx = simplify(subs(conx, y, y_sol));
                        conz = subs(conz, y, y_sol);
                        conx = subs(conx, y, y_sol);
                    end
                    conx=simplify(conx);
                    %cony=simplify(cony);
                    %conz=simplify(conz);
                    X_sol_set.(clmntocmp)=...
                            {x_sol; y_sol; z_sol; conx; cony; conz};
               %end
        
                if isempty(y_sol)
                    x_sol = [];
                    z_sol = [];
                % repeat the x and z for the different y solutions
                elseif length(y_sol) == 2
                    if length(x_sol) < 2
                        x_sol = repmat(x_sol, 2, 1);
                    end
                    if length(conx) < 2
                        conx = repmat(conx, 2, 1);
                    end
                    if length(z_sol) < 2
                        z_sol = repmat(z_sol, 2, 1);
                    end
                    if length(conz) < 2
                        conz = repmat(conz, 2, 1);
                    end
                end
            else
                if phi == pi/2
                    % solve for x
                    [x_sol, ~, conx] = solve(eq2, x,...
                            ReturnConditions=true);
                    % solve for y (maybe) in terms of z
                    [y_sol, ~, cony] = solve(...
                        subs(eq3, x, x_sol), y,...
                            ReturnConditions=true);
                    % solve for z
                    [z_sol, ~, conz] = solve(...
                        subs(eq1, {x, y}, {x_sol, y_sol}),...
                        z, ReturnConditions=true);
    
                    % no z in x, so no need for substitution of z into
                    % x
    
                    y_sol = subs(y_sol, z, z_sol);
                    cony = subs(cony, z, z_sol);
    
                    if isempty(z_sol)
                        x_sol = [];
                        y_sol = [];
                    % repeat the x and y for the possibly different z
                    % solutions
                    elseif length(z_sol) == 2
                        if length(x_sol) < 2
                            x_sol = repmat(x_sol, 2, 1);
                        end
                        if length(conx) < 2
                            conx = repmat(conx, 2, 1);
                        end
                        if length(y_sol) < 2
                            y_sol = repmat(y_sol, 2, 1);
                        end
                        if length(cony) < 2
                            cony = repmat(cony, 2, 1);
                        end
                    end
                else
                    % solve for z (maybe) in terms of x
                    [z_sol, ~, conz] = solve(eq2, z,...
                            ReturnConditions=true);
                    % solve for y (maybe) in terms of x
                    [y_sol, ~, cony] = solve(...
                        subs(eq3, z, z_sol), y,...
                            ReturnConditions=true);
                    % solve for x
                    [x_sol, ~, conx] = solve(subs(eq1, {z, y},...
                        {z_sol, y_sol}), x, ReturnConditions=true);
                    z_sol = subs(z_sol, x, x_sol);
                    y_sol = subs(y_sol, x, x_sol);
                    conz = subs(conz, x, x_sol);
                    cony = subs(cony, x, x_sol);
    
                    if isempty(x_sol)
                        z_sol = [];
                        y_sol = [];
                    % repeat the x and y for the possibly different z
                    % solutions
                    elseif length(x_sol) == 2
                        if length(z_sol) < 2
                            z_sol = repmat(z_sol, 2, 1);
                        end
                        if length(conz) < 2
                            conz = repmat(conz, 2, 1);
                        end
                        if length(y_sol) < 2
                            y_sol = repmat(y_sol, 2, 1);
                        end
                        if length(cony) < 2
                            cony = repmat(cony, 2, 1);
                        end
                    end
                end
                %conx=simplify(conx);
                cony=simplify(cony);
                conz=simplify(conz);
                X_sol_set.(clmntocmp)=...
                    {x_sol; y_sol; z_sol; conx; cony; conz};
            end
        end
        if ~isempty(x_sol)
            try
                conx = simplify(subs(conx,...
                    {Rss ds g hs o ph ps},...
                    {R_s d gamma h omega phi psi}));
                cony = simplify(subs(cony,...
                    {Rss ds g hs o ph ps},...
                    {R_s d gamma h omega phi psi}));
                conz = simplify(subs(conz,...
                    {Rss ds g hs o ph ps},...
                    {R_s d gamma h omega phi psi}));
                x_sol = x_sol(logical(conx) & logical(cony) & logical(conz));
                y_sol = y_sol(logical(conx) & logical(cony) & logical(conz));
                z_sol = z_sol(logical(conx) & logical(cony) & logical(conz));
                x_sol = subs(x_sol,...
                    {g Rss ds ph ps o},...
                    {gamma R_s d phi psi omega});
                y_sol = subs(y_sol,...
                    {g Rss ds ph ps o},...
                    {gamma R_s d phi psi omega});
                z_sol = subs(z_sol,...
                    {g Rss ds ph ps o},...
                    {gamma R_s d phi psi omega});
                x_sol = double(x_sol);
                y_sol = double(y_sol);
                z_sol = double(z_sol);
                if length(x_sol) == 2
                    % if the solutions are roughly equal, assume no
                    % intersections since the integration across the
                    % intersections is 0
                    if abs(x_sol(2) - x_sol(1)) < 1e-14 &&...
                           abs(y_sol(2) - y_sol(1)) < 1e-14 &&...
                           abs(z_sol(2) - z_sol(1)) < 1e-14
                        x_sol = [];
                        y_sol = [];
                        z_sol = [];
                    end
                end
            catch ME
                if strcmp(ME.identifier,...
                        'symbolic:kernel:DivisionByZero')
                    % if there is a division by 0, assume there is no
                    % solution
                    x_sol = [];
                    y_sol = [];
                    z_sol = [];
                end
            end
        end
    end
end

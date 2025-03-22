function [F, varargout] = ViewFac(R_s, d, phi, omega, gamma, varargin)
    % Finds the view factor for albedo between a planet and a surface
    % using what was found in 
    % https://www.sciencedirect.com/science/article/pii/S001793102033413X?fr=RR-2&ref=pdf_download&rr=8f4ba51f696f9123#fig0003
    % Hoewever, there were some errors in the article and have been
    % fixed in this program.
    % properly adjusted.
    % required inputs:
    %   R_s: Radius of full sphere that the cap is part of 
    %   d: distance of the infinitesimal sphere to the center of the
    %   sphere the cap is a part of
    %   phi: central angle between center of spherical cap to -z direction  
    %   omega: angle between the z axis and the normal vector
    %   gamma: angle between the +x direction the normal surface
    % parametric inputs:
    %   normal: infinitesimal surface normal
    %   psi: central angle between center of spherical cap to its edge
    %   X_sol_set: table of intersection solutions for the plane of the
    %   plate, projection plane, and sphere
    % outputs:
    %   F: view factor
    %   varargout: table of intersection solutions for the plane of the
    %   plate, projection plane, and sphere
    % note: psi should specified if light can't be assumed to only cover
    % half of a planet
    % note: infinitesimal surface is on the z-axis in the negative
    % direction relative to the sphere
    
    p = inputParser;
    %{
    val_R_s = @(i) valid_R_s(i);
    val_d = @(i) valid_d(i, R_s);
    val_phi = @(i) valid_zero_180(i);
    val_omega = @(i) valid_zero_180(i);
    val_gamma = @(i) valid_gamma(i);
    val_psi = @(i) valid_psi(i);
    %}
    addRequired(p, 'R_s')
    addRequired(p, 'd')
    addRequired(p, 'phi')
    addRequired(p, 'omega')
    addRequired(p, 'gamma')
    addParameter(p, 'psi', 90);
    addParameter(p, 'debug', false);
    addParameter(p, 'X_sol_set',...
        table('RowNames',{'x_sol', 'y_sol', 'z_sol', 'conx', 'cony',...
        'conz'}), @(T) isempty(T) || istable(T));
    
    parse(p, R_s, d, phi, omega, gamma, varargin{:});
    
    R_s = p.Results.R_s;
    d = p.Results.d;
    omega = p.Results.omega;
    psi = p.Results.psi;
    phi = p.Results.phi;
    gamma = p.Results.gamma;
    debug = p.Results.debug;
    X_sol_set = p.Results.X_sol_set;

    % forces the empty table to have a table with the correct rows
    if isempty(X_sol_set)
        X_sol_set = table('RowNames',...
            {'x_sol', 'y_sol', 'z_sol', 'conx', 'cony', 'conz'});
    elseif size(X_sol_set, 1) ~= 6
        X_sol_set = table('RowNames',...
            {'x_sol', 'y_sol', 'z_sol', 'conx', 'cony', 'conz'});
    end

    L3 = 0;
    
    % angle between the z axis and the tangential line from
    % the infinitesimal surface to the spherical surface
    % direction projected to the xy plane
    theta = asind(R_s./d);

    % radius of projection plane
    R = R_s.*cosd(theta);

    % correct definition for h
    % direct distance of plate to projection plane
    h = R.*cotd(theta);
    
    nom = [sind(omega).*cosd(gamma); sind(omega).*sind(gamma); cosd(omega)];
    l1 = nom(1, :);
    m1 = nom(2, :);
    n1 = nom(3, :);
    
    %{
    if or(isnan(gamma), not(isnumeric(gamma)))
        cosre = l1./sind(omega);
        sinre = m1./sind(omega);
        if sind(omega) == 0
            gamma = 0;
        elseif and(cosre >= 0, sinre >= 0)
            gamma = acosd(cosre);
        elseif and(cosre < 0, sinre >= 0)
            gamma = acosd(cosre);
        elseif and(sinre >=0, sinre < 0)
            gamma = asind(sinre);
        else
            gamma = -acosd(cosre);
        end
    end
    %}
    
    % initiate new state variables
    beta_0 = [];
    beta_1 = [];
    alpha_1 = [];
    alpha_2 = [];
    alpha_3 = [];
    X_start = [];
    X_end = [];

    % calculates view factor for each case
    if phi - psi >= 90 - theta
        F = 0; % case 1 (dont achieve for ir)
    elseif omega >= theta + 90
        F = 0; % case 2
    elseif psi - phi + theta >= 90
        if omega <= 90 - theta
            % case 3
            F = F3();
        else
            % case 4
            F = F4();
        end
    elseif phi + psi + theta <= 90
        % case 16 - 21 (dont achieve for ir or albedo typically)
        F = caseSixtTwent();
    else
        if omega <= 90 - theta
            % case 5 (dont achieve for ir)
            F = F5();
        else
            % case 6 - 15 (dont achieve for ir)
            [X_in, alphas, num_inter] = intersect_calc();
            if num_inter == 0
                % case 6 - 9 (dont achieve for ir)
                F = caseSixNine();
            elseif num_inter == 1
                % case 10 - 11 (dont achieve for ir)
                X_start = X_in.*h./X_in(3);
                alpha_2 = alphas;
                if and(gamma > 0, gamma < 180)
                    % case 10
                    F = F10();
                else
                    % case 11
                    F = F11();
                end
            elseif num_inter == 2
                % case 12 - 15 (dont achieve for ir)
                X_2_in = X_in(:, 1);
                X_3_in = X_in(:, 2);
                X_start = X_3_in.*h./X_3_in(3);
                X_end = X_2_in.*h./X_2_in(3);
                alpha_2 = alphas(1);
                alpha_3 = alphas(2);
                F = caseTwelvFift();
            end
        end
    end
    
    varargout{1} = X_sol_set;


    function [X_inter, alphas, num_inter] = intersect_calc()
        
        % gets the interserction if it exists and updates the X_sol_set
        % table
        [x_sol, y_sol, z_sol] =...
            X_intersect2(d, R_s, gamma, phi, omega, psi);
        
        X_inter = [];
        num_inter = [];
        if isempty(x_sol)
            num_inter = 0;
        elseif or(or(length(x_sol) == 2,...
                length(y_sol) == 2),...
                length(z_sol) == 2)
            num_inter = 2;
            X_inter = [x_sol.'; y_sol.'; z_sol.'];
        elseif or(or(isscalar(x_sol),...
                isscalar(y_sol)),...
                isscalar(z_sol))
            num_inter = 1;
            X_inter = [x_sol; y_sol; z_sol];
        end
        

        % debugging error to check to make sure there isn't a situation
        % with too many solutions
        if isempty(num_inter)
            error('solution has too many solutions somehow, need to fix')
        end

        alphas = [];
        if not(isempty(x_sol))
            % calculation of alphas for the intersections
            %alpha = [];
            %syms alpha real
            %assumeAlso(alpha > -180)
            %assumeAlso(alpha < 180)
            alphas = zeros([1 size(X_inter, 2)]);
            % iterates across the intersections to find the alphas
            for i = 1:size(X_inter, 2)
                cos_re = (X_inter(1,i)+R_s.*cosd(psi).*sind(phi))./...
                    (R_s.*sind(psi).*cosd(phi));
                sin_re = X_inter(2,i)./(R_s.*sind(psi));
                if cosd(phi) ~= 0
                    if cos_re - 1 > 0
                        % fixes cos_re to be 1 if it is out of bounds
                        cos_re = 1;
                    elseif cos_re + 1 < 0
                        % fixes cos_re to be -1 if it is out of bounds
                        cos_re = -1;
                    end

                    if and(cos_re >= 0, sin_re > 0)
                        alpha = acosd(cos_re);
                    elseif and(cos_re < 0, sin_re >= 0)
                        alpha = acosd(cos_re);
                    elseif and(cos_re >=0, sin_re < 0)
                        alpha = asind(sin_re);
                    else
                        alpha = -acosd(cos_re);
                    end
                else
                    alpha = asind(sin_re);
                end
                alphas(i) = alpha;
            end
            %assume(alpha, 'clear')
            if length(alphas) == 2
                if alphas(2) > alphas(1)
                    % alpha_2 > alpha_3 maps to
                    % alphas_sol(2) > alphas_sol(1)

                    % rearrange to match [alpha_2, alpha_3]
                    alphas = fliplr(alphas);
                    X_inter = fliplr(X_inter);
                end
            end
        end
    end
    %{
    function [x_sol, y_sol, z_sol] = X_intersect()
        % X_intersect finds if and where there is an intersection and creates
        % table enteries of generic solutions with the conditions for
        % intersection to exist

        if omega == 0 || omega == 180 ||...
                (gamma==0 && phi==omega &&...
                abs(d.*cosd(phi)-R_s.*cosd(psi))>1e-16)...
                || ((gamma==180 || gamma==-180) && phi==180-omega &&...
                abs(d.*cosd(phi)-R_s.*cosd(psi))>1e-16)
            % guaranteed situations to have no solution
            x_sol = [];
            y_sol = [];
            z_sol = [];
            conx = [];
            cony = [];
            conz = [];
        else
            clmntocmp = "";
            % determines gamma, sets up the column name with gamma
            if gamma == 0
                clmntocmp = clmntocmp + "g=0&";
        
                % special case can occur at phi=omega
                if phi == omega
                    clmntocmp = clmntocmp + "ph=o&";
                    if d.*cosd(phi) == R_s.*cosd(psi)
                        % spe=special case
                        clmntocmp = clmntocmp + "spe&";
                    end % other case is no solution
                else
                    clmntocmp = clmntocmp + "ph~=o&";
                end
            elseif gamma == 180 || gamma == -180
                clmntocmp = clmntocmp + "g=180|g=-180&";
        
                % special case can occur at phi=180 - omega
                if phi == 180 - omega
                    clmntocmp = clmntocmp + "ph=180-o&";
                    if abs(d.*cosd(phi)-R_s.*cosd(psi))<1e-16
                        % spe=special case
                        clmntocmp = clmntocmp + "spe&";
                    end % other case is no solution
                else
                    clmntocmp = clmntocmp + "ph~=180-o&";
                end
            else
                clmntocmp = clmntocmp + "g~=0&g~=180&g~=-180&";
            end
            % determines phi, sets up the column name with phi
            if phi == 0
                clmntocmp = clmntocmp + "ph=0&";
            elseif phi == 180
                clmntocmp = clmntocmp + "ph=180&";
            elseif phi == 90
                clmntocmp = clmntocmp + "ph=90&";
            elseif phi < 90
                clmntocmp = clmntocmp + "ph<90&";
            else %if phi > 90
                clmntocmp = clmntocmp + "ph>90&";
            end
            % determines omega, sets up the column name with omega
            if omega == 90
                clmntocmp = clmntocmp + "o=90";
            elseif omega < 90
                clmntocmp = clmntocmp + "o<90";
            else %if omega > 90
                clmntocmp = clmntocmp + "o>90";
            end
    
            x = [];
            y = [];
            z = [];
            syms x y z real
            
            gam = [];
            omeg = [];
            ps = [];
            ph = [];
            R_s_sym = [];
            d_sym = [];
            thet = [];
            h_sym = [];
            R_sym = [];
            syms gam omeg ps ph R_s_sym d_sym thet h_sym R_sym real
            if any(strcmp(clmntocmp,...
                X_sol_set.Properties.VariableNames))
                x_sol = X_sol_set{'x_sol', clmntocmp}{1};
                y_sol = X_sol_set{'y_sol', clmntocmp}{1};
                z_sol = X_sol_set{'z_sol', clmntocmp}{1};
                conx = X_sol_set{'conx', clmntocmp}{1};
                cony = X_sol_set{'cony', clmntocmp}{1};
                conz = X_sol_set{'conz', clmntocmp}{1};
            else
                
                assumeAlso(R_s_sym>0)
                assumeAlso(d_sym>R_s_sym)
                assumeAlso(z >= d_sym-R_s_sym)
                assumeAlso(thet==asind(R_s_sym/d_sym))
                assumeAlso(R_sym == R_s_sym.*cosd(thet))
                assumeAlso(h_sym == R_sym./tand(thet))
                assumeAlso(z <= h_sym)
                assumeAlso(ps>0 & ps<=90)
                
                eq1 = x.^2+y.^2+(z-d_sym).^2==R_s_sym.^2;
                eq2 = x.*sind(ph)+z.*cosd(ph)==...
                    d_sym.*cosd(ph)-R_s_sym.*cosd(ps);
                eq3 = x.*sind(omeg).*cosd(gam)+...
                y.*sind(omeg).*sind(gam)+...
                z.*cosd(omeg)==0;
        
                assumeAlso(omeg>0 & omeg<180)
                % determines gamma, sets up the column name with gamma
                if gamma == 0
                    assumeAlso(gam==0)    
                    % special case can occur at phi=omega
                    if phi == omega
                        assumeAlso(ph==omeg) % specific assumption
                        if d.*cosd(phi) == R_s.*cosd(psi)
                            assumeAlso(d_sym.*cosd(ph) == R_s_sym.*cosd(ps))
                        end % other case is no solution
                    else
                        % assumption to avoid special case
                        assumeAlso(ph~=180 - omeg)
                    end
                elseif gamma == 180 || gamma == -180
                    assumeAlso(gam==180 | gam==-180)
            
                    % special case can occur at phi=180 - omega
                    if phi == 180 - omega
                        assumeAlso(ph==180 - omeg) % specific assumption
                        if abs(d.*cosd(phi)-R_s.*cosd(psi))<1e-16
                            assumeAlso(d_sym.*cosd(ph) == R_s_sym.*cosd(ps))
                        end % other case is no solution
                    else
                        % assumption to avoid special case
                        assumeAlso(ph~=180 - omeg)
                    end
                else
                    assumeAlso(gam~=0 & gam~=180 & gam~=-180)
                end
                % determines phi, sets up the column name with phi
                if phi == 0
                    assumeAlso(ph == 0)
                elseif phi == 180
                    assumeAlso(ph == 180)
                elseif phi == 90
                    assumeAlso(ph == 90)
                elseif phi < 90
                    assumeAlso(ph>0 & ph<90)
                else %if phi > 90
                    assumeAlso(ph>90 & ph<180)
                end
                % determines omega, sets up the column name with omega
                if omega == 90
                    assumeAlso(omeg==90)
                elseif omega < 90
                    assumeAlso(omeg<90)
                else %if omega > 90
                    assumeAlso(omeg>90)
                end
                
                eq2 = simplify(eq2);
                eq3 = simplify(eq3);
            
                % gamma==0, gamma==180 and gamma==-180 have special case 
                % scenarios where the plates made from eq2 and eq3 can be
                % parallel
                if gamma==0 || (gamma==180 || gamma==-180)
                    % special case can occur at gamma=0, phi=omega and
                    % d.*cosd(phi)=R_s.*cosd(psi)
                    % special case can occur at gamma=180 or -180, phi=180-omega 
                    % and d.*cosd(phi)=R_s.*cosd(psi)
                    if (gamma==0 && phi==omega) || ((gamma==180 ||...
                            gamma==-180) && phi==180-omega)
                        if abs(d.*cosd(phi)-R_s.*cosd(psi))<1e-16
                            x_sol = solve(subs(eq3, z, h_sym), x);
                            [y_sol, pamy, cony] = solve(...
                                subs(eq1, {z, x}, {h_sym, x_sol}), y,...
                                ReturnConditions=true);
                            conx = symtrue;
                            conz = symtrue;
                            
                            % repeat the x and z for the different y
                            % solutions
                            x_sol = repmat(x_sol, 2, 1);
                            z_sol = [h_sym; h_sym];
                            conx = repmat(conx, 2, 1);
                            conz = repmat(conz, 2, 1);
                            
                            X_sol_set.(clmntocmp)=...
                                {x_sol; y_sol; z_sol; conx; cony; conz};
                        end
                    else
                        if phi == 90
                            % solve for x
                            [x_sol, pamx, conx] = solve(eq2, x,...
                                    ReturnConditions=true);
                            % solve for z
                            [z_sol, pamz, conz] = solve(...
                                subs(eq3, x, x_sol), z,...
                                    ReturnConditions=true);
                            % solve for y
                            [y_sol, pamy, cony] = solve(...
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
                            [z_sol, pamz, conz] = solve(eq2, z,...
                                    ReturnConditions=true);
                            % solve for x (maybe) in terms of y
                            [x_sol, pamx, conx] = solve(...
                                subs(eq3, z, z_sol), x,...
                                    ReturnConditions=true);
                            % gets z in terms of y (if possible)
                            z_sol = subs(z_sol, x, x_sol);
                            conz = subs(conz, x, x_sol);
                            % solve for y
                            [y_sol, pamy, cony] = solve(...
                                subs(eq1, {z, x}, {z_sol, x_sol}),...
                                y, ReturnConditions=true);
                            z_sol = subs(z_sol, y, y_sol);
                            x_sol = subs(x_sol, y, y_sol);
                            conz = simplify(subs(conz, y, y_sol));
                            conx = simplify(subs(conx, y, y_sol));
                        end
                        
                        X_sol_set.(clmntocmp)=...
                                {x_sol; y_sol; z_sol; conx; cony; conz};
                    end
            
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
                    if phi == 90
                        % solve for x
                        [x_sol, pamx, conx] = solve(eq2, x,...
                                ReturnConditions=true);
                        % solve for y (maybe) in terms of z
                        [y_sol, pamy, cony] = solve(...
                            subs(eq3, x, x_sol), y,...
                                ReturnConditions=true);
                        % solve for z
                        [z_sol, pamz, conz] = solve(...
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
                        [z_sol, pamz, conz] = solve(eq2, z,...
                                ReturnConditions=true);
                        % solve for y (maybe) in terms of x
                        [y_sol, pamy, cony] = solve(...
                            subs(eq3, z, z_sol), y,...
                                ReturnConditions=true);
                        % solve for x
                        [x_sol, pamx, conx] = solve(subs(eq1, {z, y},...
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
    
                    X_sol_set.(clmntocmp)=...
                        {x_sol; y_sol; z_sol; conx; cony; conz};
                end
            end
            if ~isempty(x_sol)
                try
                    conx = simplify(subs(conx,...
                        {gam R_s_sym d_sym ph ps omeg thet h_sym R_sym},...
                        {gamma R_s d phi psi omega theta h R}));
                    cony = simplify(subs(cony,...
                        {gam R_s_sym d_sym ph ps omeg thet h_sym R_sym},...
                        {gamma R_s d phi psi omega theta h R}));
                    conz = simplify(subs(conz,...
                        {gam R_s_sym d_sym ph ps omeg thet h_sym R_sym},...
                        {gamma R_s d phi psi omega theta h R}));
                    x_sol = x_sol(logical(conx) & logical(cony) & logical(conz));
                    y_sol = y_sol(logical(conx) & logical(cony) & logical(conz));
                    z_sol = z_sol(logical(conx) & logical(cony) & logical(conz));
                    x_sol = subs(x_sol,...
                        {gam R_s_sym d_sym ph ps omeg thet h_sym R_sym},...
                        {gamma R_s d phi psi omega theta h R});
                    y_sol = subs(y_sol,...
                        {gam R_s_sym d_sym ph ps omeg thet h_sym R_sym},...
                        {gamma R_s d phi psi omega theta h R});
                    z_sol = subs(z_sol,...
                        {gam R_s_sym d_sym ph ps omeg thet h_sym R_sym},...
                        {gamma R_s d phi psi omega theta h R});
                    x_sol = double(x_sol);
                    y_sol = double(y_sol);
                    z_sol = double(z_sol);
                    if length(x_sol) == 2
                        % if the solutions are roughly equal, choose one
                        if abs(x_sol(2) - x_sol(1)) < 1e-16 &&...
                               abs(y_sol(2) - y_sol(1)) < 1e-16 &&...
                               abs(z_sol(2) - z_sol(1)) < 1e-16
                            x_sol = x_sol(1);
                            y_sol = y_sol(1);
                            z_sol = z_sol(1);
                        end
                    end
                catch ME
                    if strcmp(ME.identifier,...
                            'symbolic:kernel:DivisionByZero')
                        % if there is a division by 0, assume there is no
                        % solution
                        conx = [];
                        cony = [];
                        conz = [];
                        x_sol = [];
                        y_sol = [];
                        z_sol = [];
                    end
                end 
            end
        end
    end
    %}
    function F = caseSixNine()
        % calculation of view factor for case 6-9

        % gets values
        if isempty(beta_1)
            beta_1 = acosd(-cosd(psi).*sind(phi)./cosd(theta)+...
            (sind(theta)-cosd(psi).*cosd(phi))./cosd(theta).*...
            cosd(phi)./sind(phi));
        end

        % decides the view factor
        if and(beta_1 >= 0, beta_1 <= 90)
            if or(and(gamma >= beta_1, gamma <= 180 - beta_1),...
                    and(gamma >= -180 + beta_1, gamma <= -beta_1))
                if h*tand(omega-90) >= R*cosd(abs(gamma)-beta_1)
                    % case 8
                    F = F8();
                elseif h*tand(omega-90) <= R*cosd(abs(gamma)+beta_1)
                    % case 7
                    F = F7();
                end
            elseif or(and(gamma > 180 - beta_1, gamma <= 180),...
                    and(gamma >= -180, gamma < -180 + beta_1))
                if h*tand(90-omega) >= R*cosd(beta_1)
                    % case 9
                    F = F9();
                else
                    % case 8
                    F = F8();
                end
            else
                if or(h*tand(omega-90) > R*cosd(beta_1),...
                        abs(h*tand(omega-90) - R*cosd(beta_1)) < 1e-14)
                    % case 6
                    F = 0;
                else
                    % case 7
                    F = F7();
                end
            end
        elseif and(beta_1 > 90, beta_1 <= 180)
            if or(and(gamma >= 180 - beta_1, gamma <= beta_1),...
                    and(gamma >= -beta_1, gamma <= -180 + beta_1))
                if or(h*tand(omega-90) > R*cosd(beta_1-abs(gamma)),...
                        abs(h*tand(omega-90) -...
                        R*cosd(beta_1-abs(gamma))) < 1e-14)
                    % case 6
                    F = 0;
                elseif h*tand(omega-90) <= R*cosd(-beta_1-abs(gamma))
                    % case 9
                    F = F9();
                end
            elseif or(and(gamma > beta_1, gamma <= 180),...
                    and(gamma >= -180, gamma < -beta_1))
                if h*tand(90-omega) <= R*cosd(beta_1)
                    % case 8
                    F = F8();
                else
                    % case 9
                    F = F9();
                end
            else
                if or(h*tand(omega-90) > R*cosd(beta_1),...
                        abs(h*tand(omega-90) - R*cosd(beta_1)) < 1e-14)
                    % case 6
                    F = 0;
                else
                    % case 7
                    F = F7();
                end
            end
        else
            error('beta_1 outside of range')
        end
    end
    function F = caseTwelvFift()
        % calculation of view factor for case 12-15
        
        if isempty(alpha_1)
            alpha_1 = acosd((sind(theta)-cosd(psi).*cosd(phi))./...
            (sind(psi).*sind(phi)));
        end
        if or(or(or(isempty(alpha_2), isempty(X_start)),...
                isempty(alpha_3)), isempty(X_end))
            [X_inter, alphas_in, ~] = intersect_calc();
            alpha_2 = alphas_in(1);
            alpha_3 = alphas_in(2);
            X_2 = X_inter(:, 1);
            X_3 = X_inter(:, 2);
            X_start = X_3.*h/X_3(3);
            X_end = X_2.*h/X_2(3);
        end
        alpha_bar = (alpha_2+alpha_3)/2;
        v = [-R_s.*cosd(psi).*sind(phi)+...
            R_s.*cosd(alpha_bar).*sind(psi).*cosd(phi);...
            R_s.*sind(alpha_bar).*sind(psi);...
            d-R_s.*cosd(psi).*cosd(phi)-...
            R_s.*cosd(alpha_bar).*sind(psi).*sind(phi)];
        if d*cosd(phi) >= R_s*cosd(psi)
            if nom.' * v > 0
                % case 12
                F = F12();
            else
                % case 13
                F = F13();
            end
        else
            if nom.' * v < 0
                % case 14
                F = F14();
            else
                % case 15
                F = F15();
            end
        end
    end

    function F = caseSixtTwent()
        if phi == 0
            theta_prime = atan2(R_s.*sind(theta), d-R_s.*cosd(psi));
            if 90+theta_prime <= omega
                % case 16
                F = 0;
            elseif and(omega >= 0, omega <= 90-theta_prime)
                % case 18
                F = F18();
            else
                % case 19
                F = F19();
            end
        else
            if omega == 0
                % case 17
                F = F17();
            else
                [X_inter, alphas_in, num_in] = intersect_calc();
                if or(num_in == 0, num_in == 1)
                    w = [-R_s.*cosd(psi).*sind(phi); 0;...
                        d-R_s.*cosd(psi).*cosd(phi)];
                    if nom.' * w > 0
                        % case 17
                        F = F17();
                    else
                        % case 16
                        F = F16();
                    end
                elseif num_in == 2
                    X_2 = X_inter(:, 1);
                    X_3 = X_inter(:, 2);
                    X_start = X_3.*h./X_3(3);
                    X_end = X_2.*h./X_2(3);
                    alpha_2 = alphas_in(1);
                    alpha_3 = alphas_in(2);
                    alpha_bar = (alpha_2+alpha_3)/2;
                    v = [-R_s.*cosd(psi).*sind(phi)+...
                        R_s.*cosd(alpha_bar).*sind(psi).*cosd(phi);...
                        R_s.*sind(alpha_bar).*sind(psi);...
                        d-R_s.*cosd(psi).*cosd(phi)-...
                        R_s.*cosd(alpha_bar).*sind(psi).*sind(phi)];
                    if nom.' * v > 0
                        % case 20
                        F = F20;
                    else
                        % case 21
                        F = F21;
                    end
                end
            end
        end
    end
    function F = F3()
        if debug
            figure(1)
            clf
        end
        F = L1(180, -180);
        if debug
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis equal
        end
        
    end

    function F = F4()
        if debug
            figure(2)
            clf
        end
        % gets any undefined values
        if isempty(beta_0)
            beta_0 = acosd(-h.*tand(90-omega)./R);
        end
        if isempty(X_start)
            X_start = [R.*cosd(gamma-beta_0); R.*sind(gamma-beta_0); h];
        end
        if isempty(X_end)
            X_end = [R.*cosd(gamma+beta_0); R.*sind(gamma+beta_0); h];
        end
        
        F = L1(gamma+beta_0, gamma-beta_0)+...
            L4(X_start, X_end);
        if debug
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis equal
        end
    end

    function F = F5()
        if debug
            figure(3)
            clf
        end
        % gets any undefined values
        if isempty(beta_1)
            beta_1 = acosd(-cosd(psi).*sind(phi)./cosd(theta)+...
            (sind(theta)-cosd(psi).*cosd(phi))./cosd(theta).*...
            cosd(phi)./sind(phi));
        end
        if isempty(alpha_1)
            alpha_1 = acosd((sind(theta)-cosd(psi).*cosd(phi))./...
            (sind(psi).*sind(phi)));
        end

        % view factor
        F = L1(2*180-beta_1, beta_1)+...
            L2(alpha_1, -alpha_1);
        if debug
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis equal
        end
    end
    
    function F = F7()
        F = F5()+F4()-F3();
    end

    function F = F8()
        F = F4();
    end

    function F = F9()
        F = F5();
    end

    function F = F10()
        if debug
            figure(4)
            clf
        end
        % gets any undefined values
        if isempty(beta_0)
            beta_0 = acosd(-h.*tand(90-omega)./R);
        end
        if isempty(beta_1)
            beta_1 = acosd(-cosd(psi).*sind(phi)./cosd(theta)+...
            (sind(theta)-cosd(psi).*cosd(phi))./cosd(theta).*...
            cosd(phi)./sind(phi));
        end
        if isempty(alpha_1)
            alpha_1 = acosd((sind(theta)-cosd(psi).*cosd(phi))./...
            (sind(psi).*sind(phi)));
        end
        if or(isempty(alpha_2), isempty(X_start))
            [X_2, alpha_2, ~] = intersect_calc();
            X_start = X_2.*h/X_2(3);
        end
        if isempty(X_end)
            X_end = [R.*cosd(gamma+beta_0); R.*sind(gamma+beta_0); h];
        end

        F = L1(gamma+beta_0, beta_1) +...
            L2(alpha_1, alpha_2)+L3+...
            L4(X_start, X_end);
        if debug
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end

    function F = F11()
        gamma = -gamma;
        nom = [sind(omega).*cosd(gamma);...
            sind(omega).*sind(gamma);...
            cosd(omega)];
        % -gamma inverts the y component of X_start and alpha_2
        X_start(2) = -X_start(2);
        alpha_2 = -alpha_2;
        F = F10();
    end

    function F = F12()
        if debug
            figure(5)
            clf
        end
        % gets any undefined values
        if or(or(or(isempty(alpha_2), isempty(alpha_3)),...
                isempty(X_start)), isempty(X_end))
            [Xs_in, alphas_in, ~] = intersect_calc();
            X_2 = Xs_in(:, 1);
            X_3 = Xs_in(:, 2);
            alpha_2 = alphas_in(1);
            alpha_3 = alphas_in(2);
            X_start = X_3.*h./X_3(3);
            X_end = X_2.*h./X_2(3);
        end

        F = L2(alpha_2, alpha_3) +...
            L3 + L4(X_start, X_end) + L3;
        if debug
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
        
    end

    function F = F13()
        F = F5() - F12();
    end

    function F = F14()
        % uses the already calculated X_start, X_end, alpha_2, alpha_3 for
        % two intersections to determine F12 first
        F12_eval = F12();
        % resets the X_start and X_end so F4 recalculates them for case 4
        X_start = [];
        X_end = [];
        F4_eval = F4();
        %{
        [Xs_in, alphas_in, ~] = intersect_calc();
        X_2 = Xs_in(:, 1);
        X_3 = Xs_in(:, 2);
        alpha_2 = alphas_in(1);
        alpha_3 = alphas_in(2);
        X_start = X_3.*h/X_3(3);
        X_end = X_2.*h/X_2(3);
        %}
        %X_start = [];
        %X_end = [];
        
        F = F4_eval - (F3() - F5() + F12_eval);
    end

    function F = F15()
        % uses the already calculated X_start, X_end, alpha_2, alpha_3 for
        % two intersections to determine F12 first
        F12_eval = F12();
        % resets the X_start and X_end so F4 recalculates them for case 4
        X_start = [];
        X_end = [];
        F4_eval = F4();
        %{
        [Xs_in, alphas_in, ~] = intersect_calc();
        X_2 = Xs_in(:, 1);
        X_3 = Xs_in(:, 2);
        alpha_2 = alphas_in(1);
        alpha_3 = alphas_in(2);
        X_start = X_3.*h/X_3(3);
        X_end = X_2.*h/X_2(3);
        %}
        F = F4_eval + F12_eval;
    end

    function F = F17()
        if debug
            figure(6)
            clf
        end
        delta = [];
        syms delta
        F = limit(L2(180-delta, -180+delta), delta, 0, 'right');
        if debug
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end        
    end

    function F = F18()
        theta_prime = atan2d(R_s.*sind(theta), d-R_s.*cosd(psi));
        R_s_prime = R_s.*sind(psi)./cosd(theta_prime);
        d_prime = d-R_s.*cosd(psi)+R_s_prime.*sind(theta_prime);
        R_prime = R_s_prime.*cosd(theta_prime);
        h_prime = R_prime./tand(theta_prime);

        theta = theta_prime;
        R_s = R_s_prime;
        d = d_prime;
        R = R_prime;
        h = h_prime;

        F = F3();
    end

    function F = F19()

        theta_prime = atan2d(R_s.*sind(theta), d-R_s.*cosd(psi));
        R_s_prime = R_s.*sind(psi)./cosd(theta_prime);
        d_prime = d-R_s.*cosd(psi)+R_s_prime.*sind(theta_prime);
        R_prime = R_s.*cosd(theta_prime);
        h_prime = R_prime./tand(theta_prime);
        beta_0_prime = acosd(-h_prime.*tand(90-omega)./R_prime);

        theta = theta_prime;
        R_s = R_s_prime;
        d = d_prime;
        R = R_prime;
        h = h_prime;
        beta_0 = beta_0_prime;
        X_start = [R_prime.*cosd(gamma-beta_0_prime);...
            R_prime.*sind(gamma-beta_0_prime); h_prime];
        X_end = [R_prime.*cosd(gamma+beta_0_prime);...
            R_prime.*sind(gamma+beta_0_prime); h_prime];

        F = F4();
    end

    function F = F20()
        if debug
            figure(7)
            clf
        end
        % gets any undefined values
        if or(or(or(isempty(alpha_2), isempty(alpha_3)),...
                isempty(X_start)), isempty(X_end))
            [Xs_in, alphas_in, ~] = intersect_calc();
            X_2 = Xs_in(:, 1);
            X_3 = Xs_in(:, 2);
            alpha_2 = alphas_in(1);
            alpha_3 = alphas_in(2);
            X_start = X_3.*h/X_3(3);
            X_end = X_2.*h/X_2(3);
        end

        F = L2(alpha_2, alpha_3) + L3 +...
            L4(X_start, X_end) + L3;
        if debug
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end        
    end


    function F = F21()
        if debug
            figure(8)
            clf
        end
        % gets any undefined values
        if or(or(or(isempty(alpha_2), isempty(alpha_3)),...
                isempty(X_start)), isempty(X_end))
            [Xs_in, alphas_in, ~] = intersect_calc();
            X_2 = Xs_in(:, 1);
            X_3 = Xs_in(:, 2);
            alpha_2 = alphas_in(1);
            alpha_3 = alphas_in(2);
            X_start = X_3.*h./X_3(3);
            X_end = X_2.*h./X_2(3);
        end

        syms delta
        F = limit(L2(180-delta, -180+delta), delta, 0,...
            'right') + L2(alpha_2, alpha_3) + L3 +...
            L4(X_start, X_end) + L3;
        if debug
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end

    function res = L1(beta_st, beta_en, varargin)
        % L1 computes the line integration of the boundary along the 
        % projection plane for the sphere in view
        % required inputs:
        %   beta_st: starting angle of boundary in radians
        %   beta_en: end angle of boundary in radians
        % inputs needed, but can maybe be defined with other inputs:
        %   nom: normal vector of the plate
        %   h: direct distance of plate to projection plane
        %   R: radius of projection plane
        % other inputs:
        %   R_s: Radius of full sphere that the cap is part of 
        %   psi: central angle between center of spherical cap to its edge
        %   phi: central angle between center of spherical cap to -z 
        %   direction
        %   theta: angle between the z axis and the tangential line
        %   from the infinitesimal surface to the spherical surface
        %   gamma: angle between the +x direction and the normal surface
        %   direction projected to the xy plane
        %   d: distance of the infinitesimal sphere to the center of the
        %   sphere the cap is a part of
        % outputs:
        %   res: the result of the line integration on the projection plane
        lp = inputParser;
        addRequired(lp, 'beta_st');
        addRequired(lp, 'beta_en');
        addParameter(lp, 'd', d);
        addParameter(lp, 'R_s', R_s);
        addParameter(lp, 'phi', phi);
        addParameter(lp, 'omega', omega);
        addParameter(lp, 'gamma', gamma);
        addParameter(lp, 'psi', psi);
        addParameter(lp, 'R', R);
        addParameter(lp, 'h', h);
        addParameter(lp, 'normal', nom);
        addParameter(lp, 'theta', theta);
        parse(lp, beta_st, beta_en, varargin{:})
        noml = lp.Results.normal;
        l1l = noml(1);
        m1l = noml(2);
        n1l = noml(3);
        Rl = lp.Results.R;
        hl = lp.Results.h;
        if debug
            hold on
            plot3(R.*cosd(linspace(beta_st, beta_en)),...
                R.*sind(linspace(beta_st, beta_en)),...
                h.*ones([1, 100]), 'r');
            text(R.*cosd(beta_st), R.*sind(beta_st), h, 'start', 'Color', 'r')
            %text(R.*cosd(beta_en), R.*sind(beta_en), h, 'end', 'Color', 'r')
            hold off
        end
        
        res = l1l.*hl.*Rl.*...
            (sind(beta_en)-sind(beta_st))./(2*pi*(Rl.^2+hl.^2))...
            +m1l.*hl.*Rl.*...
            (-cosd(beta_en)+cosd(beta_st))./(2*pi*(Rl.^2+hl.^2))...
            +n1l.*Rl.^2.*...
            (-beta_en*pi/180+beta_st*pi/180)./(2*pi*(Rl.^2+hl.^2));
    end
    
    function res = L2(alpha_st, alpha_en, varargin)
        % L2 computes the line integration along the edge of spherical cap
        % required inputs:
        %   alpha_st: starting angle of boundary in radians
        %   alpha_en: end angle of boundary in radians
        % inputs needed, but can maybe be defined with other inputs:
        %   nom: normal vector of the plate
        %   psi: central angle between center of spherical cap to its edge
        %   phi: central angle between center of spherical cap to -z
        %   direction
        %   theta: angle between the z axis and the tangential line
        %   from the infinitesimal surface to the spherical surface
        % other inputs:
        %   gamma: angle between the +x direction and the normal surface
        %   direction projected to the xy plane
        %   h: direct distance of plate to projection plane
        %   R_s: Radius of full sphere that the cap is part of
        %   R: radius of projection plane
        %   d: distance of the infinitesimal sphere to the center of the
        %   sphere the cap is a part of
        % outputs:
        %   res: the result of the line integration on the spherical cap
        lp = inputParser;
        addRequired(lp, 'alpha_st');
        addRequired(lp, 'alpha_en');
        addParameter(lp, 'd', d);
        addParameter(lp, 'R_s', R_s);
        addParameter(lp, 'phi', phi);
        addParameter(lp, 'omega', omega);
        addParameter(lp, 'gamma', gamma);
        addParameter(lp, 'psi', psi);
        addParameter(lp, 'R', R);
        addParameter(lp, 'h', h);
        addParameter(lp, 'normal', nom);
        addParameter(lp, 'theta', theta);
        parse(lp, alpha_st, alpha_en, varargin{:})
        phil = lp.Results.phi;
        psil = lp.Results.psi;
        noml = lp.Results.normal;
        l1l = noml(1);
        m1l = noml(2);
        n1l = noml(3);
        thetal = lp.Results.theta;
%{
        if all(isnan(alpha_st))
            error('Starting alpha not defined')
        end
        if all(isnan(alpha_en))
            error('End alpha not defined')
        end
    
        % parameter error checking to see if inputs are correctly defined
        if all(isnumeric(alpha_st) == 1)
            if not(all(abs(alpha_st) <= 180))
                error('Starting alpha out of bounds from -pi to pi.')
            end
        end
        if all(isnumeric(alpha_en) == 1)
            if not(all(abs(alpha_en) <= 180))
                error('Starting alpha out of bounds from -pi to pi.')
            end
        end
%}
        A = 2*sind(thetal).*sind(psil).*sind(phil)./...
            (1+sind(thetal).^2-2*sind(thetal).*cosd(psil).*cosd(phil));
        l2 = (1-sind(thetal).*cosd(psil).*cosd(phil)).*sind(thetal).*sind(psil);
        l3 = sind(thetal).^2.*sind(psil).^2.*sind(phil);
        m2 = sind(thetal).*sind(psil).*cosd(phil)-...
            sind(thetal).^2.*cosd(psil).*sind(psil);
        n2 = sind(thetal).^2.*cosd(psil).*sind(psil).*sind(phil);
        n3 = sind(thetal).^2.*sind(psil).^2.*cosd(phil);
        
        
        integ1 = -l2./A.*(alpha_en*pi/180-alpha_st*pi/180)+...
            (l2./A-l3).*(2./(1-A.^2).^(1/2)).*...
            (atan(((1+A)./(1-A)).^(1/2).*tand(alpha_en/2))-...
            atan(((1+A)./(1-A)).^(1/2).*tand(alpha_st/2)));
        
        integ2 = m2./A.*(log(1-A.*cosd(alpha_en))-...
            log(1-A.*cosd(alpha_st)));
        
        integ3 = -n2./A.*(alpha_en*pi/180-alpha_st*pi/180)+...
            (n2./A-n3).*(2./(1-A.^2).^(1/2)).*...
            (atan(((1+A)./(1-A)).^(1/2).*tand(alpha_en/2))-...
            atan(((1+A)./(1-A)).^(1/2).*tand(alpha_st/2)));
        
        
        
        res = l1l.*A./(4*pi*sind(thetal).*sind(psil).*sind(phil)).*integ1+...
            m1l.*A./(4*pi*sind(thetal).*sind(psil).*sind(phil)).*integ2+...
            n1l.*A./(4*pi*sind(thetal).*sind(psil).*sind(phil)).*integ3;

        if debug
            %integral(@(a) (l2*cosd(a)-l3)./(1-A*cosd(a)), alpha_st, alpha_en)
            %integral(@(a) (m2*sind(a))./(1-A*cosd(a)), alpha_st, alpha_en)
            %integral(@(a) (n2*cosd(a)-n3)./(1-A*cosd(a)), alpha_st, alpha_en)
            alpha = [];
            syms alpha
            x = @(alpha) -R_s.*cosd(psi).*sind(phi)+R_s.*cosd(alpha).*sind(psi).*cosd(phi);
            y = @(alpha) R_s.*sind(alpha).*sind(psi);
            z = @(alpha) d - R_s*cosd(psi)*cosd(phi)-R_s*cosd(alpha)*sind(psi)*sind(phi);
            Dx = @(a) double(subs(diff(x(alpha), alpha), alpha, a));
            Dy = @(a) double(subs(diff(y(alpha), alpha), alpha, a));
            Dz = @(a) double(subs(diff(z(alpha), alpha), alpha, a));
            int_c1 = integral(@(a) (y(a).*Dz(a)-z(a).*Dy(a))./(x(a).^2+y(a).^2+z(a).^2), alpha_st, alpha_en);
            int_c2 = integral(@(a) (x(a).*Dz(a)-z(a).*Dx(a))./(x(a).^2+y(a).^2+z(a).^2), alpha_st, alpha_en);
            int_c3 = integral(@(a) (y(a).*Dx(a)-x(a).*Dy(a))./(x(a).^2+y(a).^2+z(a).^2), alpha_st, alpha_en);
            res_c = l1./(2.*pi).*int_c1+m1./(2.*pi).*int_c2+n1./(2.*pi).*int_c3;
            %disp(abs(res_c - res) < 1e-5);
            hold on
            plot3(-R_s.*cosd(psi).*sind(phi)+R_s.*cosd(linspace(alpha_st, alpha_en))*sind(psi)*cosd(phi),...
                R_s*sind(linspace(alpha_st, alpha_en))*sind(psi),... %h*ones([1 100]), 'g')
                d - R_s*cosd(psi)*cosd(phi)-R_s*cosd(linspace(alpha_st, alpha_en))*sind(psi)*sind(phi),...
                'g')
            text(-R_s.*cosd(psi).*sind(phi)+R_s.*cosd(alpha_st)*sind(psi)*cosd(phi),...
                R_s*sind(alpha_st)*sind(psi),...
                d - R_s*cosd(psi)*cosd(phi)-R_s*cosd(alpha_st)*sind(psi)*sind(phi),...
                'start', 'Color', 'g')
            %text(-R_s.*cosd(psi).*sind(phi)+R_s.*cosd(alpha_en)*sind(psi)*cosd(phi),...
            %    R_s*sind(alpha_en)*sind(psi),...
            %    d - R_s*cosd(psi)*cosd(phi)-R_s*cosd(alpha_en)*sind(psi)*sind(phi),...
            %    'end', 'Color', 'g')
            hold off
            %res = res_c;
        end
    end

    function res = L4(X_st, X_en, varargin)
        % L4 computes the integration along a straight line on the
        % projection plane indicating the edge of the sphere out of view
        % required inputs:
        %   X_st: starting position of the sphere out of view along the
        %   projection plane
        %   X_en: end position of the boundary in the sphere out of view
        %   along the projection plane
        % inputs needed, but can maybe be defined with other inputs:
        %   nom: normal vector of the plate
        %   R_s: Radius of full sphere that the cap is part of 
        %   psi: central angle between center of spherical cap to its edge
        %   phi: central angle between center of spherical cap to -z
        %   direction
        %   theta: angle between the z axis and the tangential line
        %   from the infinitesimal surface to the spherical surface
        % other inputs:
        %   gamma: angle between the +x direction and the normal surface
        %   direction projected to the xy plane
        %   h: direct distance of plate to projection plane
        %   R: radius of projection plane
        %   d: distance of the infinitesimal sphere to the center of the
        %   sphere the cap is a part of
        % outputs:
        %   res: the result of the line integration on the edge for what is
        %   out of view along the projection plane
        
        lp = inputParser;
        addRequired(lp, 'X_st');
        addRequired(lp, 'X_en');
        addParameter(lp, 'd', d);
        addParameter(lp, 'R_s', R_s);
        addParameter(lp, 'phi', phi);
        addParameter(lp, 'omega', omega);
        addParameter(lp, 'gamma', gamma);
        addParameter(lp, 'psi', psi);
        addParameter(lp, 'R', R);
        addParameter(lp, 'h', h);
        addParameter(lp, 'normal', nom);
        addParameter(lp, 'theta', theta);
        parse(lp, X_st, X_en, varargin{:})
        noml = lp.Results.normal;
        l1l = noml(1);
        m1l = noml(2);
        n1l = noml(3);
        hl = lp.Results.h;

        x_st = X_st(1);
        y_st = X_st(2);
        x_en = X_en(1);
        y_en = X_en(2);
        % z_st and z_en are just h since the projection plane is parallel
        % to the xy plane

        x_delta = x_en - x_st;
        y_delta = y_en - y_st;

        A = x_delta.^2+y_delta.^2;
        B = x_st.*x_delta+y_st.*y_delta;
        C = x_st.^2+y_st.^2+hl.^2;
        if debug
            hold on
            plot3(x_st + x_delta.*linspace(0, 1),...
                y_st + y_delta.*linspace(0, 1),...
                X_st(3) + (X_en(3)-X_st(3)).*linspace(0, 1), 'b');
            text(x_st, y_st, X_st(3), 'start', 'Color', 'b')
            %text(x_st + x_delta, y_st + y_delta, h, 'end', 'Color', 'b')
            hold off
        end
        res = (l1l.*hl.*y_delta+...
            m1l.*-hl.*x_delta+...
            n1l.*(y_st.*x_delta-x_st.*y_delta))./...
            (2*pi*(A.*C-B.^2).^(1/2)).*...
            (atan2(A+B,(A.*C-B.^2).^(1/2))-...
            atan2(B,(A.*C-B.^2).^(1/2)));
        
    end
end

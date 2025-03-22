function [F_IR, F_A] = PlanetViewFactors(surf_pamss, planet_pams, varargin)
    % PlanetViewFactors finds the view factors on a face
    % variables:
    %   surf_pamss: struct with parameters of surface positions in sun
    %   reference frame and other surface properties
    %   planet_pams: struct of planet position in sun reference 
    %   frame and other planet properties
    % parameter variables:
    %   oribt_pams: struct of orbit parameters
    %   orbit_poss: struct with array of satellite positions and the index
    %   of the current position
    %   quats: array of thermal quaternion for each surface
    %   mRod_pams: array of modified rodrigues parameters for each surface
    %   therm_frames: array pre-defined thermal reference frames for each 
    %   surface
    
    ip = inputParser;
    addParameter(ip, 'orbit_pams', []);
    addParameter(ip, 'orbit_poss', []);
    addParameter(ip, 'quats', []);
    addParameter(ip, 'mRod_pams', []);
    addParameter(ip, 'therm_frames', []);
    parse(ip, varargin{:});
    orbit_pams = ip.Results.orbit_pams;
    orbit_poss = ip.Results.orbit_poss;
    quats = ip.Results.quats;
    mRod_pams = ip.Results.mRod_pams;
    therm_frames = ip.Results.therm_frames;
    
    % planet position and the planet to sun vector
    SP = planet_pams.pos.';
    PS = -SP;

    % surface vertices and faces
    sfaces = surf_pamss.Faces;
    svertices = surf_pamss.Vertices;

    % find normal vectors for the surface faces
    if not(isempty(surf_pamss.UserData) == 1)
        sN = surf_pamss.UserData.FaceNormals.';
    else
        sN = calcNormals(sfaces, svertices).';
    end

    % determines the centroids for the surfaces in the sun frame
    STs = surf_pamss.UserData.Centroids.';
    TSs = -STs;
    % planet to surface vectors
    PTs = PS + STs;
    TPs = -PTs;
    




    % TODO: IMPLEMENT BELOW WHEN ORBIT_PAMS, ORBIT_POSS, QUATERNIONS, AND 
    % MODIFIED RODRIGUES PARAMETERS ARE GIVEN




    % can define the nadir as the z axis

    % Can define the relative path of motion around the planet 
    % projected onto a vector perpendicular to the z axis to be the x axis.
    % If the path of motion is unknown, then use the sun as a reference to
    % project a vector onto the z axis. However, there will be a problem
    % when the plate, planet, and sun align because this vector can't be
    % computed. Consider using different method of keeping track of
    % frame.
    
    % can define the cross product between the z and x as the y axis
    
    % gamma is the angle that needs to be found with the coordinate systems
    if ~isempty(therm_frames)
        x_new = therm_frames.x;
        y_new = therm_frames.y;
        z_new = therm_frames.z;
    elseif ~isempty(mRod_pams)
        error('mRod_pams not implemented yet')
    elseif ~isempty(quats)
        error('quats not implemented yet')
    elseif ~isempty(orbit_poss)
        error('orbit_poss not implemented yet')
    elseif ~isempty(orbit_pams)
        error('orbit_pams not implemented yet')
    else
        %{
        warning("No parameter given to PlanetViewFactors. It is " + ...
            "suggested to provide one due to some singularities that " + ...
            "must be handled differently when the sun to planet " + ...
            "vector aligns with the sun to surface vector.")
        %}
        

        % find the basis axes represented in the sun reference frame

        % new z axis is the nadir from each surface
        z_new = TPs./(dot(TPs, TPs, 1).^(1/2));
        
        TSs_norm = TSs./(dot(TSs, TSs, 1).^(1/2)); % normalized surface to sun vector

        % the new x and y axis for each surface
        x_new = zeros(size(z_new));
        y_new = zeros(size(z_new));

        zv = TPs(:, all(z_new == TSs_norm, 1) | all(z_new == -TSs_norm, 1));

        if ~isempty(zv)
            % sets the new y axis to be a perpdicular projection of the 
            % the sun reference frame -z axis since it is not likely to be
            % aligned with the sun to planet vector
            % sets the new x axis to be the cross product of the new y and z
            % axis
            y_ref = repmat([0; 0; 1], 1, size(zv, 2));

            x_ref = cross(...
                y_ref(:, all(z_new == TSs_norm, 1) | all(z_new == -TSs_norm, 1)),...
                z_new(:, all(z_new == TSs_norm, 1) | all(z_new == -TSs_norm, 1)), 1);
            x_new =  x_ref ./ dot(x_ref, x_ref, 1).^(1/2);

            y_new(:, all(z_new == TSs_norm, 1) | all(z_new == -TSs_norm, 1)) = cross(...
                z_new(:, all(z_new == TSs_norm, 1) | all(z_new == -TSs_norm, 1)),...
                x_new(:, all(z_new == TSs_norm, 1) | all(z_new == -TSs_norm, 1)), 1);
        end
        
        
        zv = TPs(:, any(z_new ~= TSs_norm, 1) & any(z_new ~= -TSs_norm, 1));
        if ~isempty(zv)
            % set the new x axis to be the sun to surface vector projected 
            % to be perpendicular to the nadir
            % set the new y axis to be the cross product between the nadir 
            % and the surface to sun vector

            x_ref = TSs_norm(:, any(z_new ~= TSs_norm, 1) & any(z_new ~= -TSs_norm, 1));

            y_ref = cross(...
                z_new(:, any(z_new ~= TSs_norm, 1) & any(z_new ~= -TSs_norm, 1)),...
                x_ref(:, any(z_new ~= TSs_norm, 1) & any(z_new ~= -TSs_norm, 1)), 1);
            y_new =  y_ref ./ dot(y_ref, y_ref, 1).^(1/2);

            x_new(:, any(z_new ~= TSs_norm, 1) & any(z_new ~= -TSs_norm, 1)) = cross(...
                y_new(:, any(z_new ~= TSs_norm, 1) & any(z_new ~= -TSs_norm, 1)),...
                z_new(:, any(z_new ~= TSs_norm, 1) & any(z_new ~= -TSs_norm, 1)), 1);
        end
    end
    tole = 1e-14;
    x_new_m = dot(x_new, x_new).^(1/2);
    y_new_m = dot(y_new, y_new).^(1/2);
    z_new_m = dot(z_new, z_new).^(1/2);
    if any(abs(x_new_m - 1) > tole) || any(abs(y_new_m - 1) > tole) ||...
            any(abs(z_new_m - 1) > tole)
        error("basis vector magnitudes are too large")
    end
    x_new = x_new ./ x_new_m;
    y_new = y_new ./ y_new_m;
    z_new = z_new ./ z_new_m;
    
    % calculate parameters for the view factor
    
    % angle between the planet to surface vector and planet to sun vector
    phis = acosd(dot(PTs, repmat(PS, 1, size(PTs, 2)), 1)./...
        (dot(PTs, PTs, 1).^(1/2).*dot(PS, PS, 1).^(1/2)));

    % basis matrices
    S = eye(3);
    %{
    H = [reshape(x_new, [], 1),...
        reshape(y_new, [], 1),...
        reshape(z_new, [], 1)];
    %}
    H = zeros(length(reshape(x_new, [], 1)), 3);
    H(1:3:end, :) = x_new.';
    H(2:3:end, :) = y_new.';
    H(3:3:end, :) = z_new.';
    HS = H*S; % direction cosine matricies (for all surfaces)

    % get the representation of the surface normal vectors in the 
    % H (thermal) frame
    temp = HS*sN;
    sN_new = sN;
    sN_new(1, :) = diag(temp(1:3:end, :));
    sN_new(2, :) = diag(temp(2:3:end, :));
    sN_new(3, :) = diag(temp(3:3:end, :));
    l1s = sN_new(1, :);
    m1s = sN_new(2, :);
    
    %{
    x_new_l = reshape(x_new, [], 1);
    y_new_l = reshape(y_new, [], 1);
    z_new_l = reshape(z_new, [], 1);

    md = zeros(length(x_new_l), 1);
    d1=zeros(length(x_new_l)-1, 1);
    d2=zeros(length(x_new_l)-2, 1);
    dn1=d1;
    dn2=d2;
    md(1:3:end)=x_new(1:3:end);
    md(2:3:end)=y_new(2:3:end);
    md(3:3:end)=z_new(3:3:end);
    d1(1:3:end)=y_new(1:3:end);
    d1(2:3:end)=z_new(2:3:end);
    %d1(3:3:end)=0;
    d2(1:3:end)=z_new(1:3:end);
    %d2(2:3:end)=0;
    %d2(3:3:end)=0;
    dn1(1:3:end)=x_new(2:3:end);
    dn1(2:3:end)=y_new(3:3:end);
    dn2(1:3:end)=x_new(3:3:end);
    H2 = diag(md)+diag(d1, 1)+diag(d2, 2)+diag(dn1, -1)+diag(dn2,-2);

    S2 = eye(length(x_new_l));
    HS21 = H2^1*S2;
    HS22 = H2'*S2;
    
    sN_l = reshape(sN, [], 1);
    sN_new_l1 =HS21*sN_l;
    sN_new_l2 =HS22*sN_l;
    sN_new1 = reshape(sN_new_l1, 3, []);
    sN_new2 = reshape(sN_new_l2, 3, []);
    l1s1 = sN_new1(1, :);
    m1s1 = sN_new1(2, :);
    l1s2 = sN_new2(1, :);
    m1s2 = sN_new2(2, :);
    %}
    % angle between z axis (nadir) and normal vector
    %nadir_new = repmat([0; 0; 1], 1, size(sN_new, 2));
    omegas = acosd(dot(z_new, sN, 1)./...
        (dot(z_new, z_new, 1).^(1/2)).*dot(sN, sN, 1).^(1/2));
    
    %{
    cosre1 = l1s1./sind(omegas);
    sinre1 = m1s1./sind(omegas);
    cosre2 = l1s2./sind(omegas);
    sinre2 = m1s2./sind(omegas);
    %}
    cosre = l1s./sind(omegas);
    sinre = m1s./sind(omegas);
    tole2 = 1e10;
    cosre(cosre > 1 & cosre < 1+tole2) = 1;
    sinre(sinre > 1 & sinre < 1+tole2) = 1;
    cosre(cosre < -1 & cosre > -1-tole2) = -1;
    sinre(sinre < -1 & sinre > -1-tole2) = -1;
    
    gammas = zeros(size(omegas));
    gammas(sind(omegas) == 0) = 0;
    gammas(sind(omegas) ~= 0 & cosre >= 0 & sinre >= 0) =...
        acosd(cosre(sind(omegas) ~= 0 & cosre >= 0 & sinre >= 0));
    gammas(sind(omegas) ~= 0 & cosre < 0 & sinre >= 0) =...
        acosd(cosre(sind(omegas) ~= 0 & cosre < 0 & sinre >= 0));
    gammas(sind(omegas) ~= 0 & cosre >=0 & sinre < 0) =...
        asind(sinre(sind(omegas) ~= 0 & cosre >=0 & sinre < 0));
    gammas(sind(omegas) ~= 0 & cosre < 0 & sinre < 0) =...
        -acosd(cosre(sind(omegas) ~= 0 & cosre < 0 & sinre < 0));

    ds = dot(PTs, PTs, 1).^(1/2);
    
    R_s = planet_pams.R;
    
    % albedo view factor
    F_A = zeros(1, length(omegas));
    F_A(:) = ViewFacs2(R_s + planet_pams.TOA, ds, phis, omegas, gammas);
    
    % ir view factor
    F_IR = zeros(1, length(omegas));
    F_IR(:) = ViewFacs2(R_s, ds, zeros(1, length(omegas)), omegas, gammas);
end
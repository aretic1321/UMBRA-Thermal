function [inc_S, inc_IR, inc_A] = IncW(surf_pamss, varargin)
    % Incfacep finds the incident energies coming in a face and view
    % factors
    % variables:
    %   surf_pamss: struct with parameters of surface positions in sun
    %   reference frame and other surface properties
    % parameter variables:
    %   sun_pams: struct of sun parameters
    %   planet_pams: struct of planet position in sun reference 
    %   frame and other planet properties
    %   F_IR: array of view factors for IR
    %   F_A: array of view factors for albedo
    % IMPORTANT NOTE: planet_pams, F_IR, and F_A should be specified
    % together
    % Note 2: Sun_Int does not need to be included if planet_pams is
    % included

    ip = inputParser;
    addParameter(ip, 'sun_pams', struct([]));
    addParameter(ip, 'planet_pams', struct([]));
    addParameter(ip, 'F_IR', []);
    addParameter(ip, 'F_A', []);
    parse(ip, varargin{:});
    
    sun_pams = ip.Results.sun_pams;
    planet_pams = ip.Results.planet_pams;
    F_A = ip.Results.F_A;
    F_IR = ip.Results.F_IR;

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


    if ~isempty(fieldnames(planet_pams))
        % determine if the incident energy for IR and albedo can be
        % determined
        if (isempty(F_IR) && ~isempty(F_A)) || (~isempty(F_IR) && isempty(F_A))
            error(['If F_IR or F_A is specified',...
                'both should be given as inputs.']);
        end

        % planet position and the planet to sun vector
        SP = planet_pams.pos.';
        PS = -SP;
        
        % planet to surface vectors
        PTs = PS + STs;
        TPs = -PTs;
        
        % planet to surface vector magnitudes
        PTs_mag = dot(PTs, PTs, 1).^(1/2);
        
        % planet to sun vector magnitude
        PS_mag = dot(PS, PS, 1).^(1/2);
    
        PS_mul = repmat(PS, 1, size(PTs, 2));

        % angle between the PS and PTs
        phis = acosd(dot(PS_mul, PTs, 1)./(PS_mag.*PTs_mag));

        % If angle between the PS vector and PT vector is
        % less than 90 deg, delta = 1. Otherwise, if 
        % ||PT||*cos(90 - angle between PS vector and PT vector) > radius 
        % of the planet, delta = 1. Else, delta = 0.

        deltas = zeros(1, size(PTs, 2));

        deltas(phis < 90) = 1;
        deltas(phis >= 90 & PTs_mag.*cosd(phis-90) >= planet_pams.R) = 1;
        %deltas(phis >= 90 & PTs_mag.*cos(phis) < planet_pams.R) = 0;


        % mean distance the planet is to the sun
        mag_PS_avg = planet_pams.R_PS_avg;
    
        % solar constant and adjusted solar cosntant
        S = planet_pams.SolarIrradiance;
        S_adj = S.*mag_PS_avg.^2./PS_mag.^2;


        % albedo factor
        alb = planet_pams.Albedo;

        % computes the incident albedo for a planet
        inc_A = alb.*S_adj.*F_A;

        % computes the average incident infared emission around the planet 
        sigma = 5.670373*10^-8; % Stephan Boltzman constant (W/m^2/K)
        PIR = sigma*planet_pams.T.^4;
        inc_IR = PIR.*F_IR;
    else
        inc_IR = zeros(1, size(STs, 2));
        inc_A = zeros(1, size(STs, 2));
        deltas = ones(1, size(STs, 2));
        
        % radiaton intensity at the sun's surface or the surface emission
        H_sun = sun_pams.surf_emis;
        
        % radius of the sun
        R_sun = sun_pams.R;
        
        SSAT = surf_pamss.UserData.Reference.';

        SSAT_mag = dot(SSAT, SSAT, 1).^(1/2);
        S_adj = H_sun*R_sun.^2./(SSAT_mag).^2;
    end
    
    % surface to sun vector magnitude
    TSs_mag = dot(TSs, TSs, 1).^(1/2);
    
    % angles between the surface normal and surface to sun to vector
    alps = acosd(dot(sN, TSs, 1)./(TSs_mag));

    % determine the incident solar energy
    inc_S = zeros(1, size(STs, 2));
    inc_S(alps<90) = S_adj.*deltas(alps<90).*cosd(alps(alps<90));
end

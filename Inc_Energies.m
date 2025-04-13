%% Setup
clc
clear
close all
u = symunit;
if ispc
    addpath('Functions\', 'Input\', 'Output\');
elseif ismac || isunix
    addpath('Functions/', 'Input/', 'Output/');
else
    error('Platform not supported!');
end

debug = false; %%%%%% NOTE: CHANGE TO TRUE AND SET BREAKPOINTS TO SEE FULL VISUALS

dohot = true; %%%%% NOTE: chooses whether or not to use the hot case or cold case Beta value for testing around Earth

%% Sun Parameters
sun_pams = struct('R', 695700*10^3, 'surf_emis', 62.94*10^6, 'GM', 132712*10^6*10^9);

%% load planet table
t_planets = readtable('Planet_Values.xlsx', ReadRowNames=true, VariableNamingRule='modify');
for i = 1:length(t_planets.Properties.VariableNames)
    varname = t_planets.Properties.VariableNames{i};
    viend1 = find(varname=='_', 1, 'first');
    if ~isempty(viend1)
        t_planets.Properties.VariableNames{i} = varname(1:(viend1-1));
    end
end
t_planets.Properties.DimensionNames{1} = 'PlanetNames';

%% Spacecraft Geometry and Surface Properties

alpha_MLI = 0.1;
alpha_radiator = 0.15;
epsilon_MLI = 0.02;
epsilon_radiator = 0.8;

gtype = 'cylinder';
x_length = 2; % m (needed for cube)
z_length = 7; % m (not needed for cube) (height)
y_length = 2; % m (only needed for rectangular prisim)

if strcmp(gtype, 'cylinder')
    %%%%% NOTE: CHANGE DIVISIONS BELOW AS NEEDED %%%%%
    numsides = 8;  % Numberer of side faces

    radius = x_length/2;

    theta = linspace(0, 360, numsides+1)'; % Angular divisions
    theta(end) = [];

    % Vertices
    vertices = zeros(2*length(theta), 3);
    vertices(1:end/2,:) = [radius*cosd(theta), radius*sind(theta), repmat(-z_length*0.5, length(theta), 1)];
    vertices(end/2+1:end,:) = [radius*cosd(theta), radius*sind(theta), repmat(z_length*0.5, length(theta), 1)];
    
    
    % Bottom faces
    bottopFaces = [
        1:numsides;
        (numsides)+1:2*(numsides)];
    
    % Side faces
    sideFaces = zeros(numsides, 4);
    side_inds = 1:numsides; % side face indicies
    for i = side_inds
        sideFaces(i,:) = [i, i+1, numsides+i+1, numsides+i];
    end
    sideFaces(numsides,:) = [numsides, 1, numsides+1, 2*numsides];
    sideFaces = [sideFaces, NaN(size(sideFaces, 1), size(bottopFaces, 2) - size(sideFaces, 2))];
    
    face_elms = [bottopFaces(1, :); sideFaces; bottopFaces(2, :)];
    
    % first row corresponds to the radiator alphas, second row corresponds
    % to the MLI/coating alphas, third row corresponds to the coating
    % alphas, third row corresponds to the louver epsilons
    alphas = zeros(3, 2+numsides);

    % first row corresponds to the radiator epsilons, second row 
    % corresponds to the MLI/coating epsilons, third row corresponds to the
    % louver epsilons
    epsilons = zeros(3, 2+numsides);
    
    %%%%% NOTE: ASSUMING RADIATOR ON THE SIDE, %%%%%
    %%%%% CHANGE ANGLES BELOW AS NEEDED %%%%%
    % min angle where radiator can be placed at in deg
    theta_min = 270-45;
    % max angle wehre radiator can be placed at in deg (can be greater than 360)
    theta_max = 270+45;

    % Matched Properties to the faces or elements
    alphas(1, 1+side_inds(theta>=mod(theta_min, 360) &...
        theta<=mod(theta_max, 360))) =...
        alpha_radiator;
    alphas_MLI = alpha_radiator;
    epsilons(1, 1+side_inds(theta>=mod(theta_min, 360) &...
        theta<=mod(theta_max, 360))) =...
        epsilon_radiator;
    %alphas(1, 1) = alpha_radiator; % include if you want the ends to potential have a radiator
    %alphas(1, end) = alpha_radiator; % include if you want the ends to potential have a radiator
    %epsilons(1, 1) = epsilon_radiator; % include if you want the ends to potential have a radiator
    %epsilons(1, end) = epsilon_radiator; % include if you want the ends to potential have a radiator
    alphas(2, :) = alpha_MLI;
    epsilons(2, :) = epsilon_MLI;
    
    alphas(3, :) = alpha_radiator*0.5;
    epsilons(3, :) = epsilon_radiator*0.5;
    
elseif strcmpi(gtype, 'cube') || strcmpi(gtype, 'rectangular prisim')
    
    % Vertices of the cube
    vertices = [
        -0.5, -0.5, -0.5;  % v0
         0.5, -0.5, -0.5;  % v1
         0.5,  0.5, -0.5;  % v2
        -0.5,  0.5, -0.5;  % v3
        -0.5, -0.5,  0.5;  % v4
         0.5, -0.5,  0.5;  % v5
         0.5,  0.5,  0.5;  % v6
        -0.5,  0.5,  0.5;  % v7
    ];
    
    
    % Faces of the cube (defined by vertex indices)
    face_elms = [
        1, 2, 3, 4;  % Bottom face
        5, 6, 7, 8;  % Top face
        1, 2, 6, 5;  % Front face
        2, 3, 7, 6;  % Right face
        3, 4, 8, 7;  % Back face
        4, 1, 5, 8;  % Left face
    ];
    if strcmp(gtype, 'rectangular prisim')
        vertices(:, 1) = x_length*vertices(:, 1);
        vertices(:, 2) = y_length*vertices(:, 2);
        vertices(:, 3) = z_length*vertices(:, 3);
    else
        vertices(:, 1) = x_length*vertices(:, 1);
        vertices(:, 2) = x_length*vertices(:, 2);
        vertices(:, 3) = x_length*vertices(:, 3);
    end
    
    % first row corresponds to the radiator alphas, second row corresponds
    % to the MLI/coating alphas, third row corresponds to the coating
    % alphas, third row corresponds to the louver epsilons
    alphas = zeros(3, 6);

    % first row corresponds to the radiator epsilons, second row 
    % corresponds to the MLI/coating epsilons, third row corresponds to the
    % louver epsilons
    epsilons = zeros(3, 6);
    
    % Matched Properties to the faces or elements
    alphas(1, 3) = alpha_radiator;
    epsilons(1, 3) = epsilon_radiator;
    alphas(2, :) = alpha_MLI;
    epsilons(2, :) = epsilon_MLI;
end


if length(alphas) ~= length(face_elms) || length(epsilons) ~= length(face_elms)
    error(['Change alphas and epsilons to match the faces. ' ...
        'Use the geometry in the figure to determine where they should' ...
        'go. For a nadir pointed geometry, the x-axis faces the nadir.'])
end

colors = zeros([1 length(alphas) 3]);

figure
colors(1, :, :) = create_colors(alphas, epsilons);
patc = patch("Faces", face_elms, "Vertices", vertices,...
    'FaceColor', 'flat', 'CData', colors, 'FaceAlpha', 1);
%pat = patch("Faces", F, "Vertices", P, "FaceAlpha", 0.5);
patc = updatepatchdata(patc, 'Alphas', alphas, 'Epsilons', epsilons);
quiv = plotfacenormal(patc);
% makes the patch into a struct in case the figure gets deleted
pat = struct("Faces", patc.Faces, "Vertices", patc.Vertices, "UserData",...
    patc.UserData);

% Set axis properties for better visualization
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
view(3);  % 3D view


%% Planet Parameters
% Venus, Earth, Jupiter, Uranus
planets = ["Venus", "Uranus"];
planet = '';
p_i = 1;
for planet = planets
    if isempty(planet)
        disp('No input, try again.')
    elseif strcmpi(planet, 'Venus') == 1
        planet = 'Venus';
        planet_num = 2;
    elseif strcmpi(planet, 'Earth') == 1
        planet = 'Earth';
        planet_num = 3;
    elseif strcmpi(planet, 'Jupiter') == 1
        planet = 'Jupiter';
        planet_num = 5;
    elseif strcmpi(planet, 'Uranus') == 1
        planet = 'Uranus';
        planet_num = 7;
    elseif strcmpi(planet, 'none') == 1
        planet = 'None';
        planet_num = 0;
    else
        planet = '';
        disp('Input not accepted, try again.')
    end


    if strcmpi(planet, 'none') == 0
        R_p = t_planets{{planet}, {'EquatorialRadius'}}; % equatorial radius of the planet
        R_p_vol = t_planets{{planet}, {'MeanRadius'}}; % mean volumetric radius of the planet
        a_p = t_planets{{planet}, {'SemimajorAxis'}}; % distance of the sun to the planet
        a_p = ConvertUnit(a_p*u.au, u.m);
        e_p = t_planets{{planet}, {'OrbitEccentricity'}}; % eccentricity of planet around the sun
        albedo_p = t_planets{{planet}, {'BondAlbedo'}}; % albedo of planet
        S_p = t_planets{{planet}, {'SolarIrradiance'}}; % solar constant at the planet
        T_p = t_planets{{planet}, {'BlackBodyTemperature'}}; % black body temp of planet
        mu_p = t_planets{{planet}, {'GM'}}; % gravitational constant of planet
        omega_bar = t_planets{{planet}, {'LongitudeOfPerihelion'}}; % longitude of perihelion
        Omega = t_planets{{planet}, {'LongitudeOfAscendingNode'}}; % longitude of ascending node
        lambda_0_bar = t_planets{{planet}, {'MeanLongitude'}}; % mean longitude
        H = t_planets{{planet}, {'AltitudeInput'}}; % input altitude (default)
        H = H*10^3;
        %%%%%%% NOTE: to change the altitude, change it in
        %%%%%%% Planet_Values.xlsx file

        p_red = t_planets{{planet}, {'RedColor'}}; % planet average red color
        p_blue = t_planets{{planet}, {'BlueColor'}}; % planet average blue color
        p_green = t_planets{{planet}, {'GreenColor'}}; % % planet average green color

        % top of atmosphere (km) (is 30 km for Earth, but assumed to be the same
        % for other planets with an atmosphere)
        TOA = 30*10^3;
        
        % distance from focus of an elipse to a point on the elipse
        r_mag_elp = @(f, a, e) a.*(1 - e.^2)./(1+e.*cos(f)); % equation for elipse
        
        % distance of planet to Sun
        r_p_mag = @(f) r_mag_elp(f, a_p, e_p);
        
        f = linspace(0, 2*p_i, 10000);
        f = f(1:end-1);
        R_PS_avg = mean(r_p_mag(f)); % mean distance of the planet to the sun
        
        mu_sun = sun_pams.GM; % sun gravitational constant (m^3/s^2)
        
        n = sqrt(mu_sun./a_p.^3);
        
        tau = (omega_bar*p_i/180 - lambda_0_bar*p_i/180)./n; % time of perihelion passage
        
        if planet_num == 7 % if Uranus
            sundist_type = 'max'; % choose the farthest distance for cold case
        elseif planet_num ~= 0 && planet_num ~= 7 % if any planet that isn't Uranus
            sundist_type = 'min'; % choose the closest distance for hot case
        else
            sundist_type = 'spec ang';
        end
        if isempty(sundist_type)
            disp('No input, try again.')
        elseif strcmpi(sundist_type(1), 'a') == 1
            sundist_type = 'avg';
        elseif strcmpi(sundist_type(1:2), 'mi') == 1
            sundist_type = 'min';
        elseif strcmpi(sundist_type(1:2), 'ma') == 1
            sundist_type = 'max';
        elseif strcmpi(sundist_type(1), 's') == 1
            sundist_type = 'spec ang';
        else
            sundist_type = '';
            disp('Input not accepted, try again.')
        end
        
    
        %%%%% NOTE: use Orbit_Calc.slx to make the distance
        if strcmpi(sundist_type, 'avg')
            r_p = R_PS_avg;
            f_p = acos((a_p.*(1-e_p.^2)./r_p - 1)./e_p);
            tim = calc_time(f_p, mu_sun, a_p, e_p); % time after perihelion
        elseif strcmpi(sundist_type, 'spec ang')
            f_p = [];
            while isempty(f_p)
                f_p = double(input("Give the angle you wish the planet to be around the Sun (deg): "));
                if isempty(f_p)
                    disp('No input, try again.')
                end
            end
            f_p = f_p*p_i/180;
            tim = calc_time(f_p, mu_sun, a_p, e_p); % time after perihelion
        elseif strcmpi(sundist_type, 'min')
            tim = 0; % time after perihelion
        elseif strcmpi(sundist_type, 'max')
            tim = calc_time(p_i, mu_sun, a_p, e_p); % time after perihelion
        end
    
        auto_planet_s = '';
        while  isempty(auto_planet_s)
            %auto_planet_s = input("Do you want the planet to move automatically around the sun? (y/n): ", "s");
            auto_planet_s = 'n';
            if isempty(auto_planet_s)
                disp('No input, try again.')
            elseif strcmpi(auto_planet_s(1), 'y')
                auto_planet = 1;
            elseif strcmpi(auto_planet_s(1), 'n')
                auto_planet = 0;
            else
                auto_planet_s = '';
               disp('Input not accepted, try again.')
            end
        end
        clear auto_planet_s
        
        % radius of the planet that is being used to make calculations
        rad_p = R_p_vol; %%%%%%%%% NOTE: CHANGE AS NEEDED %%%%%%%%%
        
        planet_pams = struct('pos', [], 'R_PS_avg', R_PS_avg, 'R', rad_p,...
            'TOA', TOA, 'T', T_p, 'Albedo', albedo_p,...
            'SolarIrradiance', S_p);
    else
        tim = 0;
        auto_planet = 1;
        planet_pams = struct([]);
        tau = 0;
    end
    
    
    %% Spacecraft Orbit Parameters and (maybe) Orientation    
    if planet_num == 0
        mu_p = sun_pams.GM;
        rad_p = sun_pams.R;
    end
    R_ptosat_mag = rad_p + H; % semi-major axis of satellite (assuming circular)
    
    % angle between sun and orbital plane
    if planet_num == 0 || planet_num == 7
        Beta = 0; %%%%%% Note: change as needed
    elseif planet_num == 3
        if dohot
            Beta = 90;
        else
            Beta = 0;
        end
    else
        Beta = 10; %%%%%% Note: change as needed
    end
    
    % Define (default) orbital elements of satellite
    a_sat = R_ptosat_mag;        % Semi-major axis [AU]
    e_sat = 0;        % Eccentricity [-]
    i_sat = 0;        % Inclination [deg]
    R_sat = 0;        % RAAN [deg]
    w_sat = 0;        % Argument of periapse [deg]
    f_sat = 0;        % True Anomaly [deg]
    
    % orbital period of the satellite
    T = 2*p_i*sqrt((a_sat).^3./mu_p);
    
    % Define rotation matrices from satellite orbit frame to planet frame
    M_w_v = [cosd(w_sat), -sind(w_sat), 0;...
             sind(w_sat),  cosd(w_sat), 0;...
             0,         0,        1];
    
    M_i_v = [1, 0,         0;...
             0, cosd(i_sat), -sind(i_sat);...
             0, sind(i_sat),  cosd(i_sat)];
    
    M_R_v = [cosd(R_sat), -sind(R_sat), 0;...
             sind(R_sat),  cosd(R_sat), 0;...
             0,         0,        1];
    
    %% Sim Set up
    num_elements = 600; % number of elements to get from the dynamics
    steps_act = num_elements-1; % number of steps
    
    jtim = juliandate(2000, 1, 1, 0, 0, tim-tau);
    jtim_0 = juliandate(2000, 1, 1, 0, 0, tim-tau-0.5*T);
    while ~bdIsLoaded('Orbit_Calc')
        load_system("Orbit_Calc")
    end
    
    if planet_num ~= 0
        set_param(sprintf('Orbit_Calc/Orbit Propagator\nKepler (unperturbed)'), 'centralBody', planet);
        set_param('Orbit_Calc/Planetary Ephemeris', 'target', planet);
    else
        set_param(sprintf('Orbit_Calc/Orbit Propagator\nKepler (unperturbed)'), 'centralBody', 'Sun');
        set_param('Orbit_Calc/Planetary Ephemeris', 'target', 'Sun');
    end
    
    save_system("Orbit_Calc")
    
    % change the orbital elements so that the orbit only outputs based on beta
    if planet_num~= 0
        steps = 1;
        out = sim("Orbit_Calc.slx", T);
        h_sat = out.h_sat(1, :); % initial angular momentum vector and the orbit normal
        PSAT = out.X_PSAT(1, :); % initial planet to satellite position vector
        SP = out.X_SP(1, :); % initial sun to planet position vector
        SSAT = SP + PSAT; % initial sun to satellite vector
        SP_xy = [SP(1, 1:2), 0]; % projection of the sun to planet vector on the xy plane
        w_temp = acosd(-SP_xy*PSAT.'./(norm(-SP_xy).*norm(PSAT)));
        M_w_v_temp = [cosd(w_temp), -sind(w_temp), 0;...
             sind(w_temp),  cosd(w_temp), 0;...
             0,         0,        1];
        temp_vec = (M_w_v_temp*PSAT.').';
        ang = acosd(temp_vec*-SP_xy.'./(norm(-SP_xy).*norm(SSAT)));
        if ang < 1e-12
            if w_temp > 90
                w_sat = w_temp - 90;
            else
                w_sat = w_temp + 270;
            end
        else
                w_sat = 270 - w_temp;
        end
        ang = acosd(-SP*-SP_xy.'./(norm(-SP).*norm(-SP_xy)));
        i_sat = ang + Beta;
        clear h_sat x_sat SP SP_xy w_temp M_w_v_temp temp_vec ang
    end
    
    % determines the orbit and attitude
    steps = steps_act;
    out = sim("Orbit_Calc.slx", T);
    R_ptosat = out.X_PSAT;
    vertices_rot = out.Vertices_rot;
    vertices_adj = out.Vertices_sun;
    clear steps_act
    
    %% calc numeric
    
    inc_S = zeros(size(R_ptosat,1), size(pat.Faces, 1), 1);
    inc_IR = zeros(size(R_ptosat,1), size(pat.Faces, 1), 1);
    inc_A = zeros(size(R_ptosat,1), size(pat.Faces, 1), 1);
    F_IR = zeros(size(R_ptosat,1), size(pat.Faces, 1));
    F_A = zeros(size(R_ptosat,1), size(pat.Faces, 1));
    
    fig=figure;
    figpos = get(gcf, 'Position');
    figpos(2) = round(figpos(2)/2);
    figpos(4) = round(figpos(4)*3/2);
    set(gcf, 'Position', figpos)
    view(3)
    
    if debug
        if planet_num ~= 0
            subplot(3, 1, 1)
            view(3)
            axis equal
            xlabel('X'); ylabel('Y'); zlabel('Z');
            subplot(3, 1, 2)
            view(3)
            axis equal
            xlabel('X'); ylabel('Y'); zlabel('Z');
            subplot(3, 1, 3)
            view(3)
            axis equal
            xlabel('X'); ylabel('Y'); zlabel('Z');
        else
            subplot(2, 1, 1)
            view(3)
            axis equal
            xlabel('X'); ylabel('Y'); zlabel('Z');
            subplot(2, 1, 2)
            view(3)
            axis equal
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end
    end
    
    for step = 1:size(R_ptosat,1)%round(size(R_ptosat,1)*3/4)
        vertadjust = vertices_adj(:, :, step);
        avg_sat = mean(vertadjust, 1);
        if planet_num ~= 0
            planet_pams.pos = out.X_SP(step, :);
            
            % plot graphics
            if debug % for debugging and visualization
                subplot(3, 1, 1)
                scatter3([0, planet_pams.pos(1), avg_sat(1)],...
                    [0, planet_pams.pos(2), avg_sat(2)],...
                    [0, planet_pams.pos(3), avg_sat(3)],...
                    [500, 250, 100],... % marker sizes
                    [0.9290 0.6940 0.1250;... % color given to sun
                    p_red p_blue p_green;... % color  given to planet
                    0.5 0.5 0.5], 'filled'); % color given to satellite
                
                subplot(3, 1, 2)
                scatter3([planet_pams.pos(1), avg_sat(1)],...
                    [planet_pams.pos(2), avg_sat(2)],...
                    [planet_pams.pos(3), avg_sat(3)],...
                    [250, 100],... % marker sizes
                    [p_red p_blue p_green;... % color given to planet
                    0.5 0.5 0.5], 'filled'); % color given to satellite
                
                subplot(3, 1, 3)
            end
            h1_all = patch('Faces', pat.Faces, 'Vertices', vertadjust,...
                'FaceColor', 'flat',...
                'CData', colors(1, :, :));
            h1_all = updatepatchdata(h1_all,...
                struct(...
                'Alphas', pat.UserData.Alphas, 'Epsilons',...
                pat.UserData.Epsilons,...
                'Reference', pat.UserData.Reference + planet_pams.pos + R_ptosat(step, :)));
            if debug
                quiv = plotfacenormal(h1_all);
                refreshdata
                drawnow
            end
    
            [F_IR(step, :), F_A(step, :)] = PlanetViewFactors(h1_all, planet_pams);
            [inc_S(step, :), inc_IR(step, :), inc_A(step, :)] =...
                IncW(h1_all, 'planet_pams', planet_pams,...
                'F_IR', F_IR(step, :), 'F_A', F_A(step, :));
        else
            % plot graphics
            if debug % for debugging and visualization
                subplot(2, 1, 1)
                scatter3([0, avg_sat(1)],...
                    [0, avg_sat(2)],...
                    [0, avg_sat(3)],...
                    [250, 100],... % marker sizes
                    [0.9290 0.6940 0.1250;... % color  given to sun
                    0.5 0.5 0.5], 'filled'); % color given to satellite
                axis equal
                xlabel('X'); ylabel('Y'); zlabel('Z');
                view(3)
                subplot(2, 1, 2)
                hold on
                xlabel('X'); ylabel('Y'); zlabel('Z');
                view(3)
            end
            h1_all = patch('Faces', pat.Faces, 'Vertices', vertadjust,...
                'FaceColor', 'flat',...
                'CData', colors(1, :, :));
            h1_all = updatepatchdata(h1_all,...
                struct(...
                'Alphas', pat.UserData.Alphas, 'Epsilons',...
                pat.UserData.Epsilons,...
                'Reference', pat.UserData.Reference + R_ptosat(step, :)));
            if debug
                quiv = plotfacenormal(h1_all);
                axis equal
                hold off
            end
    
            [inc_S(step, :), inc_IR(step, :), inc_A(step, :)] =...
                IncW(h1_all, 'sun_pams', sun_pams);
        end
        cla
    end
    inc_S_avgs{p_i} = mean(inc_S, 1);
    inc_IR_avgs{p_i} = mean(inc_IR, 1);
    inc_A_avgs{p_i} = mean(inc_A, 1);
    
    if strcmpi(gtype, 'cube') || strcmpi(gtype, 'rectangular prisim')
        tab = table(inc_S_avgs{p_i}.', inc_IR_avgs{p_i}.', inc_A_avgs{p_i}.', 'VariableNames', ...
            {'Inc. Solar', 'Inc. Albedo', 'Inc. IR'}, 'RowNames', ...
            {'Zenith-Y', 'Nadir-Y', 'A-Sun-P','Ram-R', 'Sun-P', 'A-Ram-R'});
        if ispc
            savename = "Output\";
        elseif ismac || isunix
            savename = "Output/";
        else
            error('Platform not supported!');
        end
        savename = strcat(savename, sprintf('%s_incWs.xls', planet));
        writetable(tab, savename,'WriteRowNames', true, 'WriteVariableNames', true);
    end
    close(fig)
    p_i = p_i + 1;
end

%% make the name of the file to save the data
% creates the path to the file
if ispc
    savename = "Output\";
elseif ismac || isunix
    savename = "Output/";
else
    error('Platform not supported!');
end
% creates the name of the file to save the data in along the path
for planet = planets
    savename = strcat(savename, planet);
end
savename = strcat(savename, "_energies");


%% save the data
save(savename, "planets", "pat",...
    "inc_S_avgs", "inc_IR_avgs", "inc_A_avgs",...
    "inc_S", "inc_IR", "inc_A");


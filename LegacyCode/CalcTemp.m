clc
clear
close all
u = symunit;
%% Sun Parameters
sun_pams = struct('R', 695700*10^3, 'surf_emis', 62.94*10^6, 'GM', 132712*10^6*10^9);

%% load planets
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

% Matched Properties to the faces
alphas = [0.1 0.1 0.1 0.8 0.1 0.1];
epsilons = [0.02 0.02 0.02 0.05 0.02 0.02];

%{
% Vertices of the cube (side length = 1, centered at origin)
[x,y,z] = meshgrid(-2:4:2);
x = x(:);
y = y(:);
z = z(:);
K = convhull(x,y,z);
nodes = [x';y';z'];
elements = K';
fig = figure; % creates a place to put the geometry
model = createpde;
g=geometryFromMesh(model,nodes,elements);

p = pdegplot(model,"FaceLabels","on");

par1 = p.Parent;
ch = par1.Children;
for i = 1:length(ch)
    if isa(ch(i), 'matlab.graphics.primitive.Patch') == 1
        F = ch(i).Faces;
        P = ch(i).Vertices;
    end
end
clf
%close fig
%figure
%}

% Vertices of the cube (side length = 1, centered at origin)
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
vertices(:, 1) = 1*vertices(:, 1);
vertices(:, 2) = 1*vertices(:, 2);
vertices(:, 3) = 1*vertices(:, 3);

% Faces of the cube (defined by vertex indices)
face_elms = [
    1, 2, 3, 4;  % Bottom face
    5, 6, 7, 8;  % Top face
    1, 2, 6, 5;  % Front face
    2, 3, 7, 6;  % Right face
    3, 4, 8, 7;  % Back face
    4, 1, 5, 8;  % Left face
];


colors = zeros([1 length(alphas) 3]);

% Trinagulation and normal calculation
DT = delaunayTriangulation(vertices);
[F, P] = freeBoundary(DT);
TR = triangulation(F,P);
norms = faceNormal(TR);

% determine how the old faces and vertices match to the new ones
oldtonew = matchfaces(face_elms, vertices, F, P);
%{
alphas_new = zeros(size(alphas,2)*2, 1);
epsilons_new = zeros(size(epsilons,2)*2, 1);
alphas_new(:) = alphas(oldtonew);
alphas = alphas_new;
epsilons_new(:) = epsilons(oldtonew);
epsilons = epsilons_new;
%}

colors(1, :, :) = [alphas', zeros(length(alphas), 1), epsilons'];
pat = patch("Faces", face_elms, "Vertices", vertices,...
    'FaceColor', 'flat', 'CData', colors);
%pat = patch("Faces", F, "Vertices", P, "FaceAlpha", 0.5);
pat = updatepatchdata(pat, 'Alphas', alphas, 'Epsilons', epsilons);
quiv = plotfacenormal(pat);
% Set axis properties for better visualization
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
view(3);  % 3D view


%% Planet Parameters
% Venus, Earth, Jupiter, Uranus
planet = '';
while isempty(planet)
    planet = input("Venus, Earth, Jupiter, or Uranus, or None? ", "s");
    if isempty(planet)
        disp('No input, try again.')
    elseif strcmpi(planet(1), 'v') == 1
        planet = 'Venus';
        planet_num = 2;
    elseif strcmpi(planet(1), 'e') == 1
        planet = 'Earth';
        planet_num = 3;
    elseif strcmpi(planet(1), 'j') == 1
        planet = 'Jupiter';
        planet_num = 5;
    elseif strcmpi(planet(1), 'u') == 1
        planet = 'Uranus';
        planet_num = 7;
    elseif strcmpi(planet, 'none') == 1
        planet = 'None';
        planet_num = 0;
    else
        planet = '';
        disp('Input not accepted, try again.')
    end
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

    % top of atmosphere (km) (is 30 km for Earth, but assumed to be the same
    % for other planets with an atmosphere)
    TOA = 30*10^3;
    
    % distance from focus of an elipse to a point on the elipse
    r_mag_elp = @(f, a, e) a.*(1 - e.^2)./(1+e.*cos(f)); % equation for elipse
    
    % distance of planet to Sun
    r_p_mag = @(f) r_mag_elp(f, a_p, e_p);
    
    f = linspace(0, 2*pi, 10000);
    f = f(1:end-1);
    R_PS_avg = mean(r_p_mag(f)); % mean distance of the planet to the sun
    
    mu_sun = sun_pams.GM; % sun gravitational constant (m^3/s^2)
    
    n = sqrt(mu_sun./a_p.^3);
    
    tau = (omega_bar*pi/180 - lambda_0_bar*pi/180)./n; % time of perihelion passage
    
    % average, specifiy angle, min, or max
    sundist_type = '';
    while isempty(sundist_type)
        sundist_type = input("How do you want to determine the distance of "+...
            "the planet to the Sun?(avg, min, max, or spec angle) ", "s");
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
    end

    %%%%% NOTE 2: instead now use Orbit_Calc.slx to make the distance
    %%%%%%%%% NOTE: make sure to change certain things as needed %%%%%%%%
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
        f_p = f_p*pi/180;
        tim = calc_time(f_p, mu_sun, a_p, e_p); % time after perihelion
        %r_p = r_p_mag(f_p);
    elseif strcmpi(sundist_type, 'min')
        tim = 0; % time after perihelion
    elseif strcmpi(sundist_type, 'max')
        tim = calc_time(pi, mu_sun, a_p, e_p); % time after perihelion
    end

    auto_planet_s = '';
    while  isempty(auto_planet_s)
        auto_planet_s = input("Do you want the planet to move automatically around the sun? (y/n): ", "s");
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
    
    %{
    % planet position vector
    p_planet = [0 r_p 0];   %%%%%%%%% NOTE: CHANGE AS NEEDED %%%%%%%%%
    %}
    
    % radius of the planet that is being used to make calculations
    rad_p = R_p_vol; %%%%%%%%% NOTE: CHANGE AS NEEDED %%%%%%%%%
    
    %{
    [sx, sy, sz] = sphere;
    sx = sx(:);
    sy = sy(:);
    sz = sz(:);
    K = convhull(sx, sy, sz);
    elements = K';
    PP = [sx(:), sy(:), sz(:)];
    PDT = triangulation(K, PP);
    shp = alphaShape(PDT.Points(:, 1), PDT.Points(:, 2), PDT.Points(:, 3));
    pk = convhull(shp.Points);
    PDT = triangulation(pk, shp.Points);
    %[PT, PV] = freeBoundary(PDT);
    PT = PDT.ConnectivityList;
    PV = PDT.Points;
    figure
    hold on
    plnpat = patch('Faces', PT, 'Vertices', PV, 'FaceAlpha', 0.5);
    plnpat = updatepatchdata(plnpat);
    plotfacenormal(plnpat);
    hold off
    view(3)
    axis equal
    [PNs, ~] = calcNormals(PT, PV, p_planet);
    planet_pams = struct('pos', p_planet, 'R_PS_avg', R_PS_avg, 'R', rad_p,...
        'TOA', TOA, 'T', T_p, 'Albedo', albedo_p,...
        'SolarIrradiance', S_p,...
        'Faces', plnpat.Faces, 'Vertices', plnpat.Vertices,...
        'FaceNormals', plnpat.UserData.FaceNormals,...
        'Areas', plnpat.UserData.Areas,...
        'Centroids', plnpat.UserData.Centroids, ...
        'Reference',  plnpat.UserData.Reference);
    %}
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
H = 500*10^3;
%Sat_mas = 1; % 1 kg is unrealistic, but the mass wont't matter for now
%Sat_Iner = eye(3); % also not realistic, but won't matter for now

if planet_num == 0
    mu_p = sun_pams.GM;
    rad_p = sun_pams.R;
end
R_ptosat_mag = rad_p + H; % semi-major axis of satellite (assuming circular)

% angle between sun and orbital plane
Beta = 0; %%%%%% Note: change as needed

% Define (default) orbital elements of satellite
a_sat = R_ptosat_mag;        % Semi-major axis [AU]
e_sat = 0;        % Eccentricity [-]
i_sat = 0;        % Inclination [deg]
R_sat = 0;        % RAAN [deg]
w_sat = 0;        % Argument of periapse [deg]
f_sat = 0;        % True Anomaly [deg]
%M_sats = linspace(0, 360, 1000);        % Mean anomalies [deg]
%M_sats = M_sats(1:(end-1));

%R_ptosat = [-R_ptosat_mag 0 0]; % vector from planet to satellite
%thetas = linspace(0, 360, 1000).'.*pi/180; % Earth Central angle of the satellite to the sun
%R_ptosat = [R_ptosat_mag*cos(thetas), zeros(size(thetas)), R_ptosat_mag*sin(thetas)]; % vector from planet to satellite
%R_ptosat = [0 R_ptosat_mag*cos(theta) R_ptosat_mag*sin(theta)]; % vector from planet to satellite
%betasprime = 0;
%betasprime = 90;

% orbital period of the satellite
T = 2*pi*sqrt((a_sat).^3./mu_p);

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
num_elements = 1000; % number of elements to get from the dynamics
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
    ang = acosd(SP*-SP_xy.'./(norm(SP).*norm(SP_xy)));
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
%% calcs



%{
h = (H+R_p)/R_p; % rad view fac
k = R_p/(H+R_p); % thesis

ind = 0;
r_sats = zeros(3, length(M_sats));

% Convert all angles to radians
i_sat = i_sat*pi/180;
R_sat = R_sat*pi/180;
w_sat = w_sat*pi/180;
M_sats = M_sats*pi/180;
% Define rotation matrices from satellite orbit frame to reference frame
M_w_v = [cos(w_sat), -sin(w_sat), 0;...
         sin(w_sat),  cos(w_sat), 0;...
         0,         0,        1];

M_i_v = [1, 0,         0;...
         0, cos(i_sat), -sin(i_sat);...
         0, sin(i_sat),  cos(i_sat)];

M_R_v = [cos(R_sat), -sin(R_sat), 0;...
         sin(R_sat),  cos(R_sat), 0;...
         0,         0,        1];
%}
%{
for M_v = M_sats
%for f_v = f_vs
    ind = ind + 1;

    % Calculate true anomaly satellite
    f_v = calc_true_anomaly(M_v,e_sat);
    
    % Calculate satellite position in the orbital frame
    r_v = a_sat*(1-e_sat^2)/(1+e_sat*cos(f_v))*[cos(f_v);...
                                          sin(f_v);...
                                          0];

    % Calculate satellite position in the Earth reference frame
    r_sats(:, ind) = M_R_v*M_i_v*M_w_v*r_v;
end
R_ptosat = r_sats.';
%}







%{
orbitnorm = cross(R_ptosat(1, :), R_ptosat(2, :));
orbitnorm = orbitnorm / norm(orbitnorm);

% computes the average infared emission from the planet 
sigma = 5.670373*10^-8; % Stephan Boltzman constant (W/m^2/K)
PIR = sigma*planet_pams.T.^4;
verts = pat.Vertices;
for t = 1:length(M_sats)
    %hold on
    %{
    R1 = [cos(thetas(t)) 0 -sin(thetas(t));...
    0 1 0;...
    sin(thetas(t)) 0 cos(thetas(t))];
    R2 = [1 0 0;...
        0 cos(thetas(t)) -sin(thetas(t));...
        0 sin(thetas(t)) cos(thetas(t))];
    R = R1;
    %R = R2
    vertadjust = (R*verts.').' +...
        planet_pams.pos + R_ptosat(t, :);
    %}

    vertadjust = (M_R_v*M_i_v*M_w_v*verts.').' + R_ptosat(t, :) +...
        planet_pams.pos;

    colors(1, :, :) = [alphas', zeros(length(alphas), 1), epsilons'];
    patadj = patch("Faces", faces, "Vertices", vertadjust,...
        'FaceColor', 'flat', 'CData', colors);
    %pat = patch("Faces", F, "Vertices", P, "FaceAlpha", 0.5);
    patadj = updatepatchdata(patadj, 'Alphas', alphas, 'Epsilons', epsilons);
    figure
    for i = 1:size(patadj.Faces,1)
        h1 = patch('Faces', 1:size(patadj.Faces(i, :), 2),...
            'Vertices', patadj.Vertices(patadj.Faces(i, :), :),...
            'FaceColor', 'flat',...
            'CData', colors(1, i, :));
        h1 = updatepatchdata(h1,...
            struct(...
            'Alphas', patadj.UserData.Alphas(i), 'Epsilons',...
            patadj.UserData.Epsilons(i),...
            'Reference', patadj.UserData.Reference,...
            'FaceNormals', patadj.UserData.FaceNormals(i, :)), ...
            'Centroids', patadj.UserData.Centroids(i, :),...
            'Areas', patadj.UserData.Areas(i, :));
        view(3)
        lambda = acos(dot(h1.UserData.FaceNormals, orbitnorm));













        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TODO: Figure out why beta is having problems
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        betasprime = acos(dot(-1*h1.UserData.Centroids, orbitnorm) ./...
            norm(-1*h1.UserData.Centroids));
        if abs(lambda + betasprime) <= pi/2
            Fs(t, i) = planet_pams.SolarIrradiance.*cos(lambda).*cos(betasprime);
        elseif abs(lambda + betasprime) >= pi/2
            Fs(t, i) = 0;
        end
        inc_S(t, i) = Fs(t, i).*planet_pams.SolarIrradiance;

        % method from rad view fac
        TP = planet_pams.pos - h1.UserData.Centroids;
        nadir = TP ./ norm(TP);
        eta = acos(dot(nadir, h1.UserData.FaceNormals)); % angle between nadir and normal
        if and(abs(eta) > acos(1/h), abs(eta) < (pi - acos(1/h)))
            x = sqrt(h.^2-1);
            y = -x.*cot(eta);
            FsP(t, i) = 1/(pi.*h.^2)*(cos(eta).*acos(y)-x.*sin(eta).*sqrt(1-y.^2))+...
                1/pi.*atan(sin(eta).*sqrt(1-y.^2)./x);
            
        elseif and(abs(eta) >= 0, abs(eta) <= acos(1/h))
            FsP(t, i) = 1.*cos(eta)./h^2;
        else
            FsP(t, i) = 0;
        end
        PS = -planet_pams.pos;
        PT = -TP;
        Fsa(t, i) = FsP(t, i).*dot(PS, PT)./(norm(PS)*norm(PT)); % from thesis
        if i == 3
            disp(eta*180/pi)
            disp(nadir)
        end
        inc_IR(t, i) = PIR.*FsP(t, i);
        inc_A(t, i) = planet_pams.SolarIrradiance.*planet_pams.Albedo...
            .*Fsa(t,i);
        clf
    end
    close
end
inc_St = sum(inc_S, 1)./t;
inc_IRt = sum(inc_IR, 1)./t;
inc_At = sum(inc_A, 1)./t;

pow_in = 0;      %%%%% CHANGE IF NEED %%%%%%%%%%%%%%

in = sum(inc_St.'.*patadj.UserData.Areas.*patadj.UserData.Alphas) +...
    sum(inc_IRt.'.*patadj.UserData.Areas.*patadj.UserData.Epsilons) +...
    sum(inc_At.'.*patadj.UserData.Areas.*patadj.UserData.Alphas) +...
    pow_in;
syms T
assume (T, 'Real')
out = sum(patadj.UserData.Areas.'.*patadj.UserData.Epsilons.*sigma.*T.^4);
double(solve(in == out))
%}
%% calc numeric

%figure(3)
[x, y, z] = meshgrid(vertices(:, 1), vertices(:, 2), vertices(:, 3));
x = x(:);
y = y(:);
z = z(:);
%k = convhull(x, y, z);
%{
oldtonew = matchfaces(pat.Faces, pat.Vertices,...
    pattri.Faces, pattri.Vertices)

%%%%%% NOTE: FIX SO THAT THE PATCH FACES, VERTICES, AND OTHER VALUES ARE
%%%%%% SORTED SUCH THAT THEY ARE IN ORDER. mETHOD: USE A FOR LOOP TO LOOP
%%%%%% OVER THE OLD SET OF FACES, USE oldtonew == INDEX, THEN CONCATENATE
%%%%%% INTO A DIFFERENT ARRAY AND PUT INTO A STRUCT OR PATCH TO BE USED IN 
%%%%%% Incfacep
%}
inc_S = zeros(size(R_ptosat,1), size(pat.Faces, 1), 1);
inc_IR = zeros(size(R_ptosat,1), size(pat.Faces, 1), 1);
inc_A = zeros(size(R_ptosat,1), size(pat.Faces, 1), 1);
F_IR = zeros(size(R_ptosat,1), size(pat.Faces, 1));
F_A = zeros(size(R_ptosat,1), size(pat.Faces, 1));

% note: figure 4 is for determing the average position of the plate
% vertices and orientation of the plate relative to the planets position
%{
figure(4)
clf
figpos = get(gcf, 'Position');
figpos(2) = round(figpos(2)/2);
figpos(4) = round(figpos(4)*3/2);
set(gcf, 'Position', figpos)
%}


figure
figpos = get(gcf, 'Position');
figpos(2) = round(figpos(2)/2);
figpos(4) = round(figpos(4)*3/2);
set(gcf, 'Position', figpos)
view(3)
hold off
debug = false; %%%%%% NOTE: CHANGE TO TRUE AND SET BREAKPOINTS TO SEE FULL VISUALS
for step = 1:size(R_ptosat,1)%round(size(R_ptosat,1)*3/4)
    %{
    for facen = 1:size(pat.Faces, 1)
        figure(4)
        verts = pat.Vertices(pat.Faces(facen, :), :);
        
        vertadjust = verts +...
            planet_pams.pos + R_ptosat(step, :);

        subplot(3, 1, 1)
        % only use avg_sat for debugging since it may not represent the
        % actual plate centroid
        avg_sat = mean(vertadjust, 1);
        scatter3([0, planet_pams.pos(1), avg_sat(1)],...
            [0, planet_pams.pos(2), avg_sat(2)],...
            [0, planet_pams.pos(3), avg_sat(3)],...
            [500, 250, 100],... % marker sizes
            [0.9290 0.6940 0.1250;... % color given to sun
            0.8500 0.3250 0.0980;... % color  given to planet
            0.5 0.5 0.5], 'filled'); % color given to satellite
        axis equal
        xlabel('X'); ylabel('Y'); zlabel('Z');
        view(3)
        subplot(3, 1, 2)
        scatter3([planet_pams.pos(1), avg_sat(1)],...
            [planet_pams.pos(2), avg_sat(2)],...
            [planet_pams.pos(3), avg_sat(3)],...
            [250, 100],... % marker sizes
            [0.8500 0.3250 0.0980;... % color  given to planet
            0.5 0.5 0.5], 'filled'); % color given to satellite
        axis equal
        xlabel('X'); ylabel('Y'); zlabel('Z');
        view(3)
        subplot(3, 1, 3)
        hold on
        xlabel('X'); ylabel('Y'); zlabel('Z');
        h1 = patch('Faces', 1:size(verts, 1), 'Vertices', vertadjust,...
            'FaceColor', 'flat',...
            'CData', colors(1, facen, :));
        h1 = updatepatchdata(h1,...
            struct(...
            'Alphas', pat.UserData.Alphas(facen), 'Epsilons',...
            pat.UserData.Epsilons(facen), 'FaceNormals',...
            pat.UserData.FaceNormals(facen, :),...
            'Reference', pat.UserData.Reference + planet_pams.pos + R_ptosat));
        view(3)
        quiv = plotfacenormal(h1);
        axis equal;
        hold off
        
        % TODO: test PlanetViewFactors
        [F_IR(step, facen), F_A(step, facen)] = PlanetViewFactors(h1, planet_pams);

        %{
        [inc_S(facen), inc_IR(facen), inc_A(facen),...
            F_IR(facen), F_A(facen)] =...
            Incfacep(h1, planet_pams);
        %}
        clf
    end
    %}
    
    %vertadjust = vertices_rot(:, : , step) +...
        %    planet_pams.pos + R_ptosat(step, :);
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
                0.8500 0.3250 0.0980;... % color  given to planet
                0.5 0.5 0.5], 'filled'); % color given to satellite
            axis equal
            xlabel('X'); ylabel('Y'); zlabel('Z');
            view(3)
            subplot(3, 1, 2)
            scatter3([planet_pams.pos(1), avg_sat(1)],...
                [planet_pams.pos(2), avg_sat(2)],...
                [planet_pams.pos(3), avg_sat(3)],...
                [250, 100],... % marker sizes
                [0.8500 0.3250 0.0980;... % color  given to planet
                0.5 0.5 0.5], 'filled'); % color given to satellite
            axis equal
            xlabel('X'); ylabel('Y'); zlabel('Z');
            view(3)
            subplot(3, 1, 3)
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
            'Reference', pat.UserData.Reference + planet_pams.pos + R_ptosat(step, :)));
        if debug
            quiv = plotfacenormal(h1_all);
            axis equal
            hold off
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
    clf
end
inc_S_avgs = mean(inc_S, 1);
inc_IR_avgs = mean(inc_IR, 1);
inc_A_avgs = mean(inc_A, 1);

figure
hold on
for i = 1:size(pat.Faces, 1)
    plot(1:size(R_ptosat,1), F_IR(:, i),'*')
    plot(1:size(R_ptosat,1), F_A(:, i),'*')
    legs{2*i - 1} = sprintf("face %i F IR", i);
    legs{2*i} = sprintf("face %i F A", i);
end
legend(legs, 'location', 'best')


%% functions
%{
function [F_IR, F_A] = Incfacep(surf_pams, planet_pams)
% Incfacep finds the incident energies coming in a face at a specific position
% variables:
%   surf_pams: struct of surface position in sun reference frame and other
%   surface properties
%   planet_pams: struct of planet position in sun reference frame and other
%   planet properties
    
    % mean distance the planet is to the sun
    mag_PS_avg = planet_pams.R_PS_avg;
    
    % surface vertices
    vertices = surf_pams.Vertices;

    % determines the centroid of a surface
    surf_pos = mean(surf_pams.Vertices, 1);
    
    %N = calcNormals(surf_pams.Faces, surf_pams.Vertices);
    % find normal vectors for the original faces

    if not(isempty(surf_pams.UserData) == 1)
        N = surf_pams.UserData.FaceNormals;
    else
        N = calcNormals(surf_pams.Faces, surf_pams.Vertices);
    end

    
    SP = planet_pams.pos;
    PS = -SP;

    % change the basis vectors such that the z axis becomes the ET
    % vector
    PT = surf_pos - planet_pams.pos;
    z_new = PT/norm(PT);
    if all(z_new == [0 1 0])
        yv = [1 0 0];
    elseif all(z_new == [0 -1 0])
        yv = [-1 0 0]; 
    else
        yv = dot([0 1 0],z_new)*[0 1 0];
    end
    yp_new = yv - dot(yv,PT)/dot(PT, PT)*PT;
    y_new = yp_new/norm(yp_new);
    x_new = cross(y_new, z_new);
    
    % change of basis stuff
    E = eye(3);
    F = [x_new.' y_new.' z_new.'];
    P = F^-1*E;
    PT_new = (P*PT.').';
    PS_new = (P*PS.').';
    
    %R = @(phi,theta) (planet_pams.R + planet_pams.TOA).*...
    %    [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
    
    mag_R = planet_pams.R + planet_pams.TOA;
    
    N_new = (P*N.').';

    H = norm(PT) - planet_pams.R; % altitude
    lambda_0 = acos((planet_pams.R+planet_pams.TOA)/...
        (planet_pams.R+H));

    %LT_new = @(phi,theta) PT_new - R(phi, theta);
    %LS_new = @(phi,theta) PS_new - R(phi,theta);
    %n_new = @(phi, theta) R(phi, theta) / (R(phi, theta)*R(phi, theta)')^(1/2);

    mag_LT_sqrd = @(phi, theta)...
        (PT_new(1)-mag_R.*sin(theta).*cos(phi)).^2+...
        (PT_new(2)-mag_R.*sin(theta).*sin(phi)).^2+...
        (PT_new(3)-mag_R.*cos(theta)).^2;
    mag_LS_sqrd = @(phi, theta) ...
        (PS_new(1)-mag_R.*sin(theta).*cos(phi)).^2+...
        (PS_new(2)-mag_R.*sin(theta).*sin(phi)).^2+...
        (PS_new(3)-mag_R.*cos(theta)).^2;
    % pi - tau = delta where tau is from spenvis source and beta2 is from
    % notes and radiation view factors (RADIATIVE VIEW FACTORS)
    beta2 = @(phi,theta) ...
        acos((N_new(1).*-1.*(PT_new(1)-mag_R.*sin(theta).*cos(phi)) +...
        N_new(2).*-1.*(PT_new(2)-mag_R.*sin(theta).*sin(phi)) +...
        N_new(3).*-1.*(PT_new(3)-mag_R.*cos(theta)))./...
        (norm(N_new).*mag_LT_sqrd(phi,theta).^(1/2)));
    theta_ang = @(phi,theta) ...
        acos(...
        (sin(theta).*cos(phi).*(PT_new(1)-mag_R.*sin(theta).*cos(phi))+...
        sin(theta).*sin(phi).*(PT_new(2)-mag_R.*sin(theta).*sin(phi))+...
        cos(theta).*(PT_new(3)-mag_R.*cos(theta)))./...
        mag_LT_sqrd(phi,theta).^(1/2));
    sigma = @(phi,theta) ...
        acos(...
        (sin(theta).*cos(phi).*(PS_new(1)-mag_R.*sin(theta).*cos(phi))+...
        sin(theta).*sin(phi).*(PS_new(2)-mag_R.*sin(theta).*sin(phi))+...
        cos(theta).*(PS_new(3)-mag_R.*cos(theta)))./...
        (mag_LS_sqrd(phi,theta).^(1/2)));
    
    % incident albedo
    F_A = integral2(@(phi,theta)...
        planet_pams.Albedo .* P_r(phi,theta).*Psi(phi,theta),...
        0, 2*pi, 0, lambda_0);

    % incident ir
    F_IR = integral2(@(phi,theta) (f_ir(phi,theta).*2.*pi),...
        0, 2*pi, 0, pi);
    

    function out = f_ir(phi,theta)
        out = zeros(size(theta, 1), size(phi, 2));
        
        cmp = cos(theta_ang(phi, theta))>0;

        out(cmp) = max(cos(beta2(phi(cmp),theta(cmp))), 0).*...
            (2.*pi.*sin(theta_ang(phi(cmp),theta(cmp)))).^-1;
    end

    function out = P_i(phi, theta)
        %cmp = theta_ang(phi,theta) < pi/2;
        %cmp = ones(size(out));
        out = max(cos(sigma(phi,theta)), 0).*mag_PS_avg.^2.*...
            mag_LS_sqrd(phi,theta).^-1;
    end
    
    function out = P_r(phi, theta)
        out = zeros(size(theta, 1), size(phi, 2));
        cmp = theta_ang(phi,theta) < pi/2;
        out(cmp) = P_i(phi(cmp), theta(cmp)) ./ (2*pi*...
            sin(theta_ang(phi(cmp), theta(cmp))));
    end

    function out = Psi(phi, theta)
        A = getarea(vertices);
        out = A.*max(cos(beta2(phi,theta)), 0) .*...
            mag_LT_sqrd(phi,theta).^-1;
    end
end
%}



function [inc_S, inc_IR, inc_A, F_IR, F_A] =...
    Incfacep(surf_pams, planet_pams)
% Incfacep finds the incident energies coming in a face and view factors
% variables:
%   surf_pams: struct of surface position in sun reference frame and other
%   surface properties
%   planet_pams: struct of planet position in sun reference frame and other
%   planet properties
    inc_S = 0;
    inc_IR = 0;
    inc_A = 0;
    F_IR = 0;
    F_A = 0;

    % computes the average infared emission from the planet 
    sigma = 5.670373*10^-8; % Stephan Boltzman constant (W/m^2/K)
    PIR = sigma*planet_pams.T.^4;
    
    % mean distance the planet is to the sun
    mag_PS_avg = planet_pams.R_PS_avg;
    
    % planet position and the planet to sun vector
    SP = planet_pams.pos
    PS = -SP
    
    % surface vertices and faces
    sfaces = surf_pams.Faces;
    svertices = surf_pams.Vertices;
    pfaces = planet_pams.Faces;
    pvertices = planet_pams.Vertices;
    % find normal vectors for the planet faces
    if not(isempty(surf_pams.UserData) == 1)
        pN = planet_pams.FaceNormals
    else
        [pN, ~] = calcNormals(pfaces, pvertices)
    end
    % find normal vectors for the surface faces
    if not(isempty(surf_pams.UserData) == 1)
        sN = surf_pams.UserData.FaceNormals
    else
        sN = calcNormals(sfaces, svertices)
    end
    
    pareas = planet_pams.Areas;
    sareas = surf_pams.UserData.Areas;
    totparea = 0;
    totsarea = 0;

    % determines the centroid for the surface
    STs = surf_pams.UserData.Centroids
    % vector from planet to surface
    PTs = STs - SP
    TSs = -STs
    % loops over the surfaces from planet_pams
    for li = 1:size(pfaces , 1)
        pverts = pvertices(planet_pams.Faces(li, :), :);
        SL = planet_pams.Centroids(li, :)

        % planet to planet surface vector
        PL = SL - SP

        LS = PS - PL
        
        lN = pN(li, :)

        cossigma = LS*lN.'./...
            (sum(LS.*LS, 2)).^(1/2)

        LTs = PTs-PL

        % determine the cos(theta) angle between the planet surface
        % normal and the vector from the planet surface to the surface
        % in space
        costhetas = LTs*lN.'./...
            (sum(LTs.*LTs, 2)).^(1/2)
        tis_tru = costhetas > 0;
        
        LTs = LTs(tis_tru, :);
        costhetas = costhetas(tis_tru);
        tareas = sareas(tis_tru);

        % loops over the surfaces from surf_pams
        if not(isempty(LTs))
            
            for ti = 1:length(find(tis_tru == 1))
                tN = sN(ti, :);
                eta = acos(-LTs(ti, :)*tN.' ./ norm(-LTs(ti, :)))
                if abs(eta) > pi/2
                    continue
                end
                
                f_IR = max(cos(eta), 0).*costhetas(ti) ./...
                    (pi.*norm(LTs(ti, :)).^2);
                f_reflect = max(cossigma, 0).*mag_PS_avg^2.*norm(LS)^-2; % *mag_PS_avg^2*norm(LS)^2?
                f_albedo = f_reflect.*f_IR;
                F_IR = F_IR + f_IR .* tareas(ti).* pareas(li);
                F_A = F_A + f_albedo .* tareas(ti).* pareas(li);
                
                inc_IR = inc_IR + PIR.*f_IR .* tareas(ti) .* pareas(li);
                inc_A = inc_A + planet_pams.SolarIrradiance .*...
                    f_albedo .* planet_pams.Albedo .*...
                    tareas(ti) .* pareas(li);
            end
            totparea = totparea + pareas(li);
        end
    end

    betas = zeros(size(surf_pams.Faces, 1), 1);
    alphas = zeros(size(surf_pams.Faces, 1), 1);
    deltas = zeros(size(surf_pams.Faces, 1), 1);
    phis = zeros(size(surf_pams.Faces, 1), 1);
    inc_Ss = zeros(size(surf_pams.Faces, 1), 1);
    % gets the parameters for calculating solar energy
    for i = 1:size(surf_pams.Faces, 1)
        TS = TSs(i, :);
        % get the angle of the sun to the normal vector of the surface
        if isfield(surf_pams.UserData, 'phi') == 0
            phis(i) = acos(dot(sN(i, :), TS)/...
                (norm(sN(i, :))*norm(TS)));
        else
            phis(i) = surf_pams.phi(i);
        end

        if not(all(planet_pams.pos == [0 0 0]))
            if isfield(surf_pams, 'UserData') == 1
                if isfield(surf_pams.UserData, 'Betas') == 0
                    betas(i) = acos(dot(TS,PS)/...
                        (norm(TS)*norm(SP)));
                else
                    betas(i) = surf_pams.UserData.Betas(i);
                end
            else
                betas(i) = acos(dot(TS,PS)/...
                    (norm(TS)*norm(SP)));
            end
            if isfield(planet_pams, 'Alphas') == 0
                alphas(i) = atand(planet_pams.R/...
                    (norm(PS)^2-planet_pams.R^2)^(1/2));
            else
                alphas(i) = planet_pams.Alphas(i);
            end
            if and(betas(i)<alphas(i),...
                    norm(TS)>(norm(PS)^2-planet_pams.R^2)^(1/2))
                deltas(i) = 0;
            else
                deltas(i) = 1;
            end
        else
            deltas(i) = 1;
        end
        % incident solar flux or illumination
        inc_Ss(i) = planet_pams.SolarIrradiance .*deltas(i).*...
            max(cos(phis(i)), 0).*...
            mag_PS_avg^2./(norm(TS).^2);
        inc_S = inc_S + inc_Ss(i).*sareas(i);
        if not(inc_Ss(i) == 0)
            totsarea = totsarea + sareas(i);
        end
    end
    if not(totsarea == 0)
        inc_S = inc_S ./ totsarea;
    else
        inc_S = 0;
    end
end


function fir_i_in_sec = matchfaces(fir_faces, fir_vertices, sec_faces, sec_vertices)
% matchfaces matches the first face indicies and puts them in an array
% where each index in the array corrosponds to the index of the second
% face.
    if size(fir_faces,2) > size(sec_faces,2)
        % identifies that a face in the first set has more vertices
        fir_big = true;
        %{
        big_verts = fir_vertices;
        big_faces = fir_faces;
        smal_verts = sec_vertices;
        smal_faces = sec_faces;
        %}
    else
        % identifies that a face in the first set has the same or less
        % amount of vertices
        fir_big = false;
        %{
        big_verts = sec_vertices;
        big_faces = sec_faces;
        smal_verts = fir_vertices;
        smal_faces = fir_faces;
        %}
    end

    if size(fir_faces,1) < size(sec_faces,1)
        less_faces = 1; % less faces in first set
    else
        less_faces = 2; % more or eqaul num of faces in first set assuming it can't have more
    end
    fir_i_in_sec = zeros(size(sec_faces,1), 0);
    js = 1:size(sec_faces,1);
    for i = 1:size(fir_faces,1) % loops across faces
        for j = js % loops across faces
        %for j = 1:size(sec_faces,1) % loops across faces
            vid1 = fir_faces(i, :);
            vid2 = sec_faces(j, :);
            
            if fir_big
                big_verts = fir_vertices(vid1, :);
                smal_verts = sec_vertices(vid2, :);
            else
                big_verts = sec_vertices(vid2, :);
                smal_verts = fir_vertices(vid1, :);
            end
            % finds if the first and second face have matching vertices or
            % if they exist in the same plane
            if all(ismember(smal_verts, big_verts, "rows") == 1)
                js = js(not(j == js));
                fir_i_in_sec(j) = i;
                break
            elseif less_faces == 1
                % gets vectors between vertices
                vec_chk = big_verts(1:size(smal_verts, 1), :) -...
                    smal_verts;

                % gets the normal vector of the small face
                [N, ~] = calcNormalSimple(smal_verts, 2);
                
                % if the dot product of the normal and vectors between
                % vertices is 0, then the faces are coplanar
                if all(vec_chk*N.' == 0)
                    js = js(not(j == js));
                    fir_i_in_sec(j) = i;
                    break
                end
            end
        end
    end
end
%{
function idx = findMatchingPoint(point, pointList)
    % Function to check if a point matches any point in a list
    % 
    % Inputs:
    % - point: 1xN array representing the point to search for
    % - pointList: MxN matrix of points to search within
    % 
    % Output:
    % - idx: Index of the matching point in pointList (0 if no match found)

    % Ensure the input is valid
    if size(pointList, 2) ~= length(point)
        error('The dimensions of the point and pointList must match.');
    end

    % Find rows in pointList that match the given point
    matches = ismember(pointList, point, 'rows');
    
    % Get the index of the match
    idx = find(matches, 1, 'first'); % Find the first match (if any)
    
    if isempty(idx)
        idx = 0; % Return 0 if no match is found
    end
end
%}
function pat = updatepatchdata(pat, varargin)
% updates patch data if a change is specified
    p = inputParser;
    addParameter(p, 'Areas', [])
    addParameter(p, 'Alphas', [])
    addParameter(p, 'Epsilons', [])
    addParameter(p, 'FaceNormals', [])
    addParameter(p, 'Reference', [])
    addParameter(p, 'Centroids', [])
    parse(p, varargin{:});
    
    areas = p.Results.Areas;
    alphas = p.Results.Alphas;
    epsilons = p.Results.Epsilons;
    norms = p.Results.FaceNormals;
    ref = p.Results.Reference;
    centroids = p.Results.Centroids;

    numfaces = size(pat.Faces, 1);

    % for each face, add associated parameters

    if isempty(pat.UserData) == 1
        defNaN = NaN(numfaces, 1);
        pat.UserData = struct('Areas', defNaN,...
            'Alphas', defNaN, 'Epsilons', defNaN, ...
            'FaceNormals', NaN(numfaces, 3), 'Reference',...
            NaN(1, 3), 'Centroids', NaN(numfaces, 3));
    end

    if not(isempty(areas))
        % sets the areas
        pat.UserData.Areas = areas;
    elseif any(isnan(pat.UserData.Areas))
        % sets the areas if they weren't given and not already set
        for i = 1:numfaces
            pat.UserData.Areas(i) =...
                getarea(pat.Vertices(pat.Faces(i, :), :));
        end
    end

    if not(isempty(alphas))
        % sets the alphas
        pat.UserData.Alphas = alphas;
    elseif any(isnan(pat.UserData.Alphas))
        % sets the alphas if they weren't given and not already set
        warning("Alphas not specified in updatepatchdata after" +...
            " initilization. Default values set to 0.")
        pat.UserData.Alphas = zeros(numfaces, 1);
    end

    if not(isempty(epsilons))
        pat.UserData.Epsilons = epsilons;
    elseif any(isnan(pat.UserData.Epsilons))
    % sets the epsilons if they weren't given and not already set
        warning("Epsilons not specified in updatepatchdata after" +...
            " initilization. Default values set to 0.")
        pat.UserData.Epsilons = zeros(numfaces, 1);
    end
    
    if not(isempty(centroids))
        % sets the centroids
        pat.UserData.Centroids = centroids;
    elseif any(isnan(pat.UserData.Centroids))
        % sets the centroids if they weren't given and not already set
        pat.UserData.Centroids = zeros(numfaces, 3);
        for i = 1:numfaces
            pat.UserData.Centroids(i,:) =...
                shapecenter(pat.Vertices(pat.Faces(i, :), :));
        end
    end

    if not(isempty(ref))
        pat.UserData.Reference = ref;
    elseif any(isnan(pat.UserData.Reference))
        % sets the reference if it wasn't given and not already set
        if numfaces == 1
            warning("Patch only has a single face. Reference not" + ...
                " specified in updatepatchdata after initilization." + ...
                " Default reference is a row vector of NaN.")
            pat.UserData.Reference = ref;
        else
            % sets the reference to the centroid of the shape
            pat.UserData.Reference = centroid(pat.Vertices);
        end
    end

    if not(isempty(norms))
        % sets the normal vectors
        pat.UserData.FaceNormals = norms;
    elseif any(isnan(pat.UserData.FaceNormals))
        % sets the normal vectors if they weren't given and not already set
        if size(pat.Faces, 1) > 1
            pat.UserData.FaceNormals =...
                calcNormals(pat.Faces, pat.Vertices, pat.UserData.Reference);
        else
            pat.UserData.FaceNormals =...
                calcNormalSimple(pat.Vertices, 2, pat.UserData.Reference);
        end
    end

    % calculates the normal vectors and sets the patch normal vectors
    if isempty(pat.FaceNormals) == 1
        parent_object = pat.Parent;
        parent_type = parent_object.Type;
        iter = 0;
        while strcmp(parent_type, 'axes') == 0  &...
                strcmp(parent_type, 'figure') == 0 & iter < 10
            parent_object = parent_object.Parent;
            parent_type = parent_object.Type;
            iter = iter + 1;
        end

        if and(not(strcmp(parent_type, 'figure') == 1),...
                not(strcmp(parent_type, 'axes') == 1))
            error('Could not reach axes object');
        end
        
        set(pat, 'FaceNormals', pat.UserData.FaceNormals)
    end
    
end

function A = getarea(vertices, varargin)
% getarea computes the area of a polygon assuming the vertices are listed
% along the rows unless otherwise specific. dim specifies the dimension
% that the coordinates of a vertex are listed across

    p = inputParser;
    addRequired(p, 'vertices')
    addOptional(p, 'dim', 2)
    parse(p, vertices, varargin{:});
    
    dim = p.Results.dim;

    if dim == 2
        v = vertices - mean(vertices, 1);
        A = area3D(v(:, 1), v(:, 2), v(:, 3));
    elseif dim == 1
        v = vertices - mean(vertices, 1);
        A = area3D(v(1, :), v(2, :), v(3, :));
    end
end

function q = plotfacenormal(p)
% input: p is a patch object
    
    cent = zeros(size(p.Faces, 1), 3);
    for facen = 1:size(p.Faces, 1)
        % get the indices of the vertices for the face
        v_is = p.Faces(facen, :);

        % compute the centroid of the shape
        cent(facen, :) = shapecenter(p.Vertices(v_is, :));
    end
    if not(isempty(p.UserData))
        norms = p.UserData.FaceNormals;
    else
        norms = calcNormals(p.Faces, p.Vertices, p.UserData.Reference);
    end
    hold on
    q = quiver3(cent(:,1), cent(:,2), cent(:,3),...
        norms(:,1), norms(:,2), norms(:,3));
    hold off
end

function [normals, order] = calcNormals(faces, vertices, varargin)
% calcNormals finds normal vectors of a face or faces. To find the normal
% vector of a single face, just make sure faces has some array where the
% size of the rows is 1.
    
    p = inputParser;
    addRequired(p, 'Faces')
    addRequired(p, 'Vertices')
    addOptional(p, 'Reference', NaN(1, 3))
    parse(p, faces, vertices, varargin{:})
    
    ref = p.Results.Reference;
    
    if size(faces, 1) == 1
        [normals, order] = calcNormalSimple(vertices, 2, ref);
    else
        DT = delaunayTriangulation(vertices);
        [T, X] = freeBoundary(DT);
        TR = triangulation(T, X);
        Ns = faceNormal(TR);
        
        % uses matchFaces to determine which normal vecotrs match to
        % which face, then applies them to the associated original face
        normals = zeros(size(faces,1), 3);
        oldtonew = matchfaces(faces,vertices, T, X);
        order = oldtonew;
        for i = 1:size(faces, 1)
            n_i = find(oldtonew == i, 1, 'first');
            normals(i, :) = Ns(n_i, :);
        end
    end
end

function [N, order] = calcNormalSimple(vertices, dim, varargin)
% calcNormalSimple finds a normal vector to a set of vertices, where the
% coordinates of each vertex is along dim (1 for row, 2 for column).
% A reference point must be included if you want the normal to face away
% from the point. If there are more than three vertices, it assumes all
% vertices are in the same plane and the shape is convex.
% Also outputs the order of the vertices that would calculate the normals
% anti-clockwise from the perspective of the reference if a reference is
% given. Otherwise, the order is the original order.
    
    p = inputParser;
    addRequired(p, 'vertices')
    addRequired(p, 'dim')
    if dim == 1
        defref = NaN(3, 1);
    elseif dim == 2
        defref = NaN(1, 3);
    end
    addOptional(p, 'Reference', defref)
    parse(p, vertices, dim, varargin{:});
    
    ref = p.Results.Reference;
    
    mid = shapecenter(vertices, dim);

    % get the vertex coordinates
    if dim == 1
        v1 = vertices(:, 1);
        v2 = vertices(:, 2);
        v3 = vertices(:, 3);
    elseif dim == 2
        v1 = vertices(1, :);
        v2 = vertices(2, :);
        v3 = vertices(3, :);
    end
    % compute two edge vectors
    e1 = v2 - v1;
    e2 = v3 - v2;
    
    % compute the normal vector
    n = cross(e1, e2);
    
    
    if dim == 1
        order = 1:size(vertices, 2);
    else
        order = 1:size(vertices, 1);
    end

    % normalize the normal vector
    if any(isnan(ref))
        N = n ./ norm(n);
    else
        ch = dot(n, mid - ref, dim); % sign gives direction
        if ch < 0
            multip = -1;
            order = fliplr(order);
        else
            multip = 1;
        end
        
        N = n ./ norm(n) .* multip;
    end
end

function mid = shapecenter(vertices, varargin)
% finds the center of a shape assuming the vertices are coplanar and the
% shape is convex
    p = inputParser;
    addRequired(p, 'vertices')
    addOptional(p, 'dim', 2)
    parse(p, vertices, varargin{:});
    dim = p.Results.dim;
    
    if dim == 1
        oppdim = 2;
    elseif dim == 2
        oppdim = 1;
    end

    % performs triangulation around a mean vertex
    if dim == 1
        vertices = vertices.';
    end
    meanvert = mean(vertices, oppdim);
    C = cat(oppdim, meanvert, vertices);
    X = zeros(size(C, 1)-1, 3);
    for i = 1:(size(C, 1)-2)
        X(i, :) = [1 i+1 i+2];
    end
    X(size(C, 1)-1, :) = [1 size(C, 1) 2];

    % computes the center of the triangles and areas
    cents = zeros(size(X, 1), 3);
    areas = zeros(size(X, 1), 1);
    for i = 1:size(X, 1)
        cents(i, :) = mean(C(X(i, :),:), oppdim);
        areas(i) = getarea(C(X(i, :), :));
    end

    % computes the center point of the shape
    mid = sum(cents.*areas, 1)./sum(areas);
    if dim == 1
        mid = mid.';
    end
end


% Function to convert mean anomaly to true anomaly
function f = calc_true_anomaly(M,e)

    % Define convergence tolerance
    tol = 1e-6;

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
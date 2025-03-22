clc;
clear;
close all;
% solar constant calculator

u = symunit;

% determine the radiation exitance at the sun's surface
S_max_earth = 1419 * u.W/u.m^2; % max solar constant from new smad
S_min_earth = 1317 * u.W/u.m^2; % min solar constant from new smad
mu_earth = 3.986004356*10^14 *u.m^3/u.s^2;
f = linspace(0, 2*pi, 10000);
f = f(1:end-1);
a_earth = 149.598*10^6* u.km;
a_earth = ConvertUnit(a_earth, u.m);
e_earth = 0.01671022;
%%%%%%%%%%%%%% change it so the radius is calculated not as an
%%%%%%%%%%%%%% approximation
r = @(f, a, e) a*(1 - e^2)./(1+e*cos(f));
r_earth = @(f) r(f, a_earth, e_earth);
mean_d_earth = mean(r_earth(f)); % mean distance of the earth to the sun
r_sun = 695700 * u.km; % mean radius of the sun
r_sun = ConvertUnit(r_sun, u.m);
S_max_sun = S_max_earth*4*pi*mean_d_earth^2/(4*pi*r_sun^2);
S_min_sun = S_min_earth*4*pi*mean_d_earth^2/(4*pi*r_sun^2);
S_max_sun = ConvertUnit(S_max_sun, u.W/u.m^2);
S_min_sun = ConvertUnit(S_min_sun, u.W/u.m^2);
Ss_sun = struct('min', S_min_sun, 'max', S_max_sun);

% determine the solar constant / irradiance max and min of a planet
a_uranus = 2867.043 * 10^6 * u.km;
a_uranus = ConvertUnit(a_uranus, u.m);
e_uranus = 0.0469;
r_uranus = @(f) r(f, a_uranus, e_uranus);
min_d_uranus = r_uranus(0);
max_d_uranus = r_uranus(pi);
S_max_uranus = S_max_sun * (4*pi*r_sun^2)/(4*pi*min_d_uranus^2);
S_min_uranus = S_min_sun * (4*pi*r_sun^2)/(4*pi*max_d_uranus^2);

alpha = 0.1;
epsilon = 0.02;
A_p = 1;
A_r = 6;
sigma = 5.670373 * 10 ^-8; % Stephan Boltzman constant (W/m^2/K)
T = @(S) (S*(alpha/epsilon)*(A_p/A_r)/sigma) ^ (1/4);
T_max_uranus = T(S_max_uranus) + 274.15;
T_min_uranus = T(S_min_uranus) + 274.15;
side1 = struct('name', 'Sun', 'A', 25, 'alpha', 0.1, 'epsilon', 0.1);
side2 = struct('name', 'A-Sun', 'A', 25, 'alpha', 0.1, 'epsilon', 0.1);
side3 = struct('name', 'side1', 'A', 25, 'alpha', 0.1, 'epsilon', 0.1);
side4 = struct('name', 'side2', 'A', 25, 'alpha', 0.1, 'epsilon', 0.1);
side5 = struct('name', 'side3', 'A', 25, 'alpha', 0.1, 'epsilon', 0.1);
side6 = struct('name', 'side4', 'A', 25, 'alpha', 0.1, 'epsilon', 0.1);

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
vertices(:, 1) = 10*vertices(:, 1);
vertices(:, 2) = 10*vertices(:, 2);
vertices(:, 3) = 10*vertices(:, 3);

DT = delaunayTriangulation(vertices);
[F, P] = freeBoundary(DT);

% Faces of the cube (defined by vertex indices)
faces = [
    1, 2, 3, 4;  % Bottom face
    5, 6, 7, 8;  % Top face
    1, 2, 6, 5;  % Front face
    2, 3, 7, 6;  % Right face
    3, 4, 8, 7;  % Back face
    4, 1, 5, 8;  % Left face
];

% Create the cube using the patch function
figure;
%alphas = linspace(0, 1, 6);
%alphas = duplicateprop(alphas);
%epsilons = linspace(1, 0, 6);
%epsilons = duplicateprop(epsilons);

alphas = [0.1 0.1 0.1 0.8 0.1 0.1];
epsilons = [0.02 0.02 0.02 0.05 0.02 0.02];
%[Fs, slist] = CheckOutwardNormals(F, P);
%Fs = fliplr(Fs);

TR = triangulation(F,P);
norms = faceNormal(TR);
oldtonew = matchfaces(faces, vertices, F, P);
alphas_new = zeros(1, size(alphas,2)*2);
epsilons_new = zeros(1, size(epsilons,2)*2);
alphas_new = alphas(oldtonew);
alphas = alphas_new;
epsilons_new = epsilons(oldtonew);
epsilons = epsilons_new;
colors = zeros([1 length(alphas) 3]);
colors(1, :, :) = [alphas', zeros(length(alphas), 1), epsilons'];
h = patch('Faces', F, 'Vertices', P, 'FaceColor', 'flat', 'CData', colors,...
    'FaceNormals', calcNormals(F, P));
h.FaceNormalsMode = 'auto';
updatepatchdata(h, alphas, epsilons)


% Set axis properties for better visualization
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
view(3);  % 3D view

%plotfacenormal(h);

H = 7800*10^3;
R_earth = 6371*10^3;
R_jupiter = 69911*10^3;
R_venus = 6051.8*10^3;
R_uranus = 25362*10^3;
R_sattoe_mag = R_uranus + H;
R_sattoe = [0 -R_sattoe_mag 0];
%scale_height = 8.5*10^3;
p_planet= [2867.043*10^6*10^3 0 0];
%low_temp = 248;
%high_temp = 306;
planet_pams = struct('pos', p_planet, 'R', R_uranus,...
    'd_atm', 0,...
    'T', 58.1, 'Albedo', 0.300);
h.set('Visible', false)
%{
[X_sph, Y_sph, Z_sph] = sphere;
hold on
surf(R_earth* X_sph+p_earth(1), R_earth* Y_sph+p_earth(2),...
    R_earth* Z_sph+p_earth(3))
hold off
%}
axis equal


for facen = 1:size(h.Faces, 1)
    verts = h.Vertices(h.Faces(facen, :), :);
    hold on
    h1 = patch('Faces', [1 2 3], 'Vertices',...
        verts +...
        planet_pams.pos + R_sattoe,...
        'FaceColor', 'flat',...
        'CData', colors(1, facen, :));
    plotfacenormal(h1);
    %axis equal
    hold off
    [inc_S(facen), inc_A(facen), inc_IR(facen)] = Incfacep(Ss_sun.max, h1)%,...
        %'planet_params', planet_pams);
end
Q_in = 0;
Q_out = 0;
syms T
for facen = 1:size(h.Faces, 1)
    Q_in = Q_in + (h.UserData.Areas(facen).*...
        (h.UserData.Alphas(facen).*(inc_A(facen)+inc_S(facen))+...
        h.UserData.Epsilons(facen).*inc_IR(facen)));
    Q_out = Q_out + h.UserData.Epsilons(facen).*sigma.*T.^4;
end
Q_env = 0;
Q_power = 500;

Tem = -1*double(solve(Q_out == Q_in + Q_env + Q_power, T, 'PrincipalValue', true))

%{
r = a*(1 - e^2)./(1+e*cos(f));


L_low = 3.828*10^26; % solar luminosity (Watts)
r_close = 1*u.au; % distance from sun in units found
r_close = ConvertUnit(r_close, u.m); % distance to sun (m)
r_far = 4583190000*u.mi; % distance from sun in units found
r_far = ConvertUnit(r_far, u.m); % distance to sun (m)



S = @(r) L_low ./ (4 * pi * r.^2);
S_close = S(r_close); % solar irradiance where uranus is at it's closest
S_far = S(r_far); % solar irradiance where uranus is at it's farthest


S_all = S(r);
S_mean = mean(S_all)
S_2mean = S(mean(r))


sigma = 5.670373 * 10 ^-8; % Stephan Boltzman constant (W/m^2/K)
R_sun = 695700 * u.km;
R_sun = ConvertUnit(R_sun, u.m);
T=(L_low/(4*pi*695700*sigma))^(1/4)
%}

function T = tempcalc1(Ss_sun, planet_orb_pams, surfs)
% variables:
%   Ss_sun: struct of min and max surface emission values of the sun
%   planet_orb_pams: struct of semimajor axis and eccentricity
%   surfs: surface parameters
    r = @(f) planet_orb_pams.a*(1 - planet_orb_pams.e^2)./...
        (1+planet_orb_pams.e*cos(f));
    r_sun = 695700000;
    r_close = r(0);
    r_far = r(pi);
    %{
    r = @(f) planet_pams.a*(1 - planet_pams.e^2)./...
            (1+planet_pams.e*cos(f));
    r_sun = 695700000;
        
    if type == 'max'
        f = 0;
        r_planet = r(f);
        S = S_sun.max * (4*pi*r_sun^2)/(4*pi*r_planet^2);
    elseif type == 'min'
        f = pi;
        r_planet = r(f);
        S = S_sun.min * (4*pi*r_sun^2)/(4*pi*r_planet^2);
    end
    %}
    S_max = Ss_sun.max * (4*pi*r_sun^2)/(4*pi*r_close^2);
    S_min = Ss_sun.min * (4*pi*r_sun^2)/(4*pi*r_far^2);
    
    T = @(S) (S*(alpha/epsilon)*(A_p/A_r)/sigma) ^ (1/4);
end

function [inc_S, inc_A, inc_IR]  = Incfacep(S_sun, surf_pams, varargin)
% Incfacep finds the incident energies coming in a face at a specific position
% variables:
%   S_sun: surface emission value of the sun
%   planet_pams: struct of planet position in sun reference frame and other
%   planet properties
%   surf_pams: struct of surface position in sun reference frame and other
%   surface properties
    
    p = inputParser;
    addParameter(p, 'planet_params', struct('pos', [0 0 0]))
    parse(p, varargin{:});
    
    planet_pams = p.Results.planet_params;
    
    r_sun = 695700000;
    surf_pos = mean(surf_pams.Vertices, 1);

    % solar constant at the distance of the surface or planet
    if all(planet_pams.pos == [0 0 0])
        S = S_sun * (4*pi*r_sun^2)/(4*pi*norm(surf_pos)^2);
    else
        S = S_sun * (4*pi*r_sun^2)/(4*pi*norm(planet_pams.pos)^2);
    end
    
    N = calcNormals(surf_pams.Faces, surf_pams.Vertices);
    
    % vector from the surface to the sun
    TS = -1*surf_pos;
    % vector from the planet to the sun
    PS = -1*planet_pams.pos;

    % get the angle of the sun to the normal vector of the surface
    if isfield(surf_pams.UserData, 'phi') == 0
        phi = acos(dot(N, TS)/...
            (norm(N)*norm(TS)));
    else
        phi = surf_pams.phi;
    end

    % gets the incident solar energy
    if not(all(planet_pams.pos == [0 0 0]))
        if isfield(surf_pams, 'UserData') == 1
            if isfield(surf_pams.UserData, 'betas') == 0
                beta = acos(dot(TS,PS)/...
                    norm(TS)*norm(planet_pams.pos));
            else
                beta = surf_pams.UserData.betas;
            end
        else
            beta = acos(dot(TS,PS)/...
                norm(TS)*norm(planet_pams.pos));
        end
        if isfield(planet_pams, 'alphas') == 0
            alpha = atand(planet_pams.R/...
                (norm(PS)^2-planet_pams.R^2)^(1/2));
        else
            alpha = surf_pams.alphas;
        end
        if beta<alpha && norm(TS)>(norm(PS)^2-planet_pams.R^2)^(1/2)
            delta = 1;
        else
            delta = 0;
        end
    else
        delta = 1;
    end

    % incident solar flux or illumination
    inc_S = S*delta*max(cos(phi), 0)*...
        norm(planet_pams.pos)^2/(norm(surf_pos)^2);


    if not(all(planet_pams.pos == [0 0 0]))
        surf_pos = mean(surf_pams.Vertices, 1);
        planet_pos = planet_pams.pos;
        PS = -planet_pos;
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
        F = [x_new' y_new' z_new'];
        P = F^-1*E;
        PT_new = (P*PT')';
        PS_new = (P*PS')';
        
        %R = @(phi,theta) (planet_pams.R + planet_pams.d_atm).*...
        %    [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
        
        mag_R = planet_pams.R + planet_pams.d_atm;
        
        N_new = (P*N')';
        
        if isfield(planet_pams, 'PIR') == 1
            PIR = planet_pams.PIR;
        else
            if isa(planet_pams.T, 'function_handle')
                PIR=@(theta_ang) 5.6703*10^-8*planet_pams.T(theta_ang).^4;
            else
                PIR = 5.6703*10^-8*planet_pams.T.^4;
            end
        end
        

        H = norm(PT) - planet_pams.R; % altitude
        lambda_0 = acos((planet_pams.R+planet_pams.d_atm)/...
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
        tau = @(phi,theta) ...
            acos((N_new(1).*(PT_new(1)-mag_R.*sin(theta).*cos(phi)) +...
            N_new(2).*(PT_new(2)-mag_R.*sin(theta).*sin(phi)) +...
            N_new(3).*(PT_new(3)-mag_R.*cos(theta)))./...
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

        % incident ir
        inc_IR = integral2(@(phi,theta) (I_E(phi,theta)),...
            0, 2*pi, 0, pi);

        % incident albedo
        inc_A = integral2(@(phi,theta)...
            planet_pams.Albedo .* P_r(phi,theta).*Psi(phi,theta),...
            0, 2*pi, 0, lambda_0);
    else
        inc_IR = 0;
        inc_A = 0;
    end

    function out = I_E(phi,theta)
        temp = zeros(size(theta, 1), size(phi, 2));
        out = zeros(size(theta, 1), size(phi, 2));
        
        cmp = cos(theta_ang(phi, theta))>0;

        temp(cmp) = max(-cos(tau(phi(cmp),theta(cmp))), 0).*...
            (2.*pi.*sin(theta_ang(phi(cmp),theta(cmp)))).^-1;
        A = zeros(size(theta, 1), size(phi, 2));
        A(cmp) = mag_R.*sin(theta(cmp));
        if isa(PIR, 'function_handle')
            out(cmp) = PIR(theta_ang(phi(cmp),theta(cmp))).*...
                temp(cmp).*A(cmp);
        else
            out = PIR.*temp.*A;
        end
    end

    function out = P_i(phi, theta)
        out = zeros(size(theta, 1), size(phi, 2));
        %cmp = theta_ang(phi,theta) < pi/2;
        %cmp = ones(size(out));
        out = S*max(cos(sigma(phi,theta)), 0)*norm(PS_new)^2.*...
            mag_LS_sqrd(phi,theta).^-1;
    end
    
    function out = P_r(phi, theta)
        out = zeros(size(theta, 1), size(phi, 2));
        cmp = theta_ang(phi,theta) < pi/2;
        out(cmp) = P_i(phi(cmp), theta(cmp)) ./ (2*pi*...
            sin(theta_ang(phi(cmp), theta(cmp)))); 
    end

    function out = Psi(phi, theta)
        out = zeros(size(theta, 1), size(phi, 2));
        A = mag_R^2.*sin(theta_ang(phi,theta));
        out = A.*max(-cos(tau(phi,theta)), 0) .*...
            mag_LT_sqrd(phi,theta).^-1;
    end
end

function p = updatepatchdata(p, alphas, epsilons)
    if isempty(p.FaceNormals) == 1
        parent_object = p.Parent;
        parent_type = parent_object.Type;
        iter = 0;
        while strcmp(parent_type, 'axes') == 0  &...
                strcmp(parent_type, 'figure') == 0 & iter < 10
            parent_object = parent_object.Parent;
            parent_type = parent_object.Type;
            iter = iter + 1;
        end

        if and(not(strcmp(parent_type, 'figure') == 1),...
                not(strcmp(parent_type, 'figure') == 1))
            error('Could not reach axes object');
        end
        set(p, 'FaceNormals', calcNormals(p.Faces, p.Vertices))
    end
    % for each face, add associated alphas, and epsilons
    p.UserData = struct('Areas', zeros(1, size(p.Faces, 1)),...
        'Alphas', alphas, 'Epsilons', epsilons);
    for i = 1:size(p.Faces, 1)
        p.UserData.Areas(i) = getarea(p.Vertices(p.Faces(i, :), :));
    end
end

function old_i_in_new = matchfaces(og_faces, og_verticies, new_faces, new_verticies)
    
    old_i_in_new = zeros(1, size(new_faces,1));
    for i = 1:size(og_faces,1)
        for j = 1:size(new_faces,1)
            count = 0;
            % check if the new face is in the old face, record it, 
            % move on if it is
            for k = 1:size(new_faces, 2)
                vidx1 = new_faces(j, k);
                
                vidx2 = og_faces(i, :);
                if not(findMatchingPoint(new_verticies(vidx1, :),...
                        og_verticies(vidx2, :)) == 0)
                    count = count + 1;
                end
                if count == 3
                    old_i_in_new(j) = i;
                end
            end
        end
    end
end
function A = getarea(vertices)
% getarea computes the area of a polygon assuming the vertices are listed
% vertically
    
    % vectors from one vertex to all other verticies
    vecs = zeros(size(vertices, 1)-1, 3);
    
    % computes the area of regular polygons
    for vn = 1:size(vecs,1)
        vecs(vn, :) = vertices(vn+1, :) - vertices(1, :);
    end
    A = 0;
    for en=1:size(vecs, 1)-1
        A = A + 1/2*norm(cross(vecs(en, :), vecs(en+1,:)));
    end
end

function prop = duplicateprop(prop)
    temp = prop;
    if size(prop, 1) == 1
        prop = zeros([1 length(prop)*2]);
    else
        prop = zeros([length(prop)*2, 1]);
    end
    prop(1:2:(end-1)) = temp;
    prop(2:2:end) = temp;
end
%{
doesn't work
function old_i_in_new = GetAssocFace(og_faces, og_verts, new_faces, new_verts)
    old_i_in_new = zeros(1, size(new_faces, 1));
    % find normal vectors for the original faces
    for facen = 1:size(og_faces, 1)
        % get the indices of the vertices for the face
        v_is = og_faces(facen, :);
        
        % get the vertex coordinates
        v1 = og_verts(v_is(1), :);
        v2 = og_verts(v_is(2), :);
        v3 = og_verts(v_is(3), :);
        
        % compute two edge vectors
        e1 = v2 - v1;
        e2 = v3 - v1;
        
        % compute the normal vector
        n = cross(e1, e2);
        
        % normalize the normal vector
        og_normals(facen, :) = n / norm(n);
    end

    % find normal vectors for the new faces
    for facen = 1:size(new_faces, 1)
        % get the indices of the vertices for the face
        v_is = new_faces(facen, :);
        
        % get the vertex coordinates
        v1 = new_verts(v_is(1), :);
        v2 = new_verts(v_is(2), :);
        v3 = new_verts(v_is(3), :);
        
        % compute two edge vectors
        e1 = v2 - v1;
        e2 = v3 - v1;
        
        % compute the normal vector
        n = cross(e1, e2);
        
        % normalize the normal vector
        new_normals(facen, :) = n / norm(n);
    end
    for newfacen = 1:size(new_faces, 1)
        for oldfacen = 1:size(og_faces, 1)
            % cross product of normal vectors
            cro = cross(og_normals(oldfacen, :), new_normals(newfacen, :));

            og_v_is = og_faces(oldfacen, :);
            new_v_is = new_faces(newfacen, :);

            % dot product of normal to the difference in verticies of the
            % original and new verticies
            pln = dot(og_normals(oldfacen, :),...
                (og_verts(og_v_is(1))-new_verts(new_v_is(1))));

            if all(cro == [0 0 0]) & pln == 0
                old_i_in_new(newfacen) = oldfacen;
            end
        end
    end
end
%}

function q = plotfacenormal(p)
% input: p is a patch object
    
    cent = zeros(size(p.Faces, 1), 3);
    for facen = 1:size(p.Faces, 1)
        % get the indices of the vertices for the face
        v_is = p.Faces(facen, :);

        % compute the centroid of the shape
        cent(facen, :) = mean(p.Vertices(v_is, :), 1);
    end

    prop = get(p);
    norms = calcNormals(p.Faces, p.Vertices);
    hold on
    q = quiver3(cent(:,1), cent(:,2), cent(:,3),...
        norms(:,1), norms(:,2), norms(:,3));
    hold off
end

function normals = calcNormals(faces, vertices)
% find normal vectors
    TR = triangulation(faces, vertices);
    normals = faceNormal(TR);
end

function out = temperVariation(low, high, sigma)
    out = zeros(size(sigma)) + low; 
    out(sigma < pi/2) = out(sigma < pi/2) +...
        (high - low) .* cos(sigma(sigma < pi/2));
end



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
    idx = find(matches, 1); % Find the first match (if any)
    
    if isempty(idx)
        idx = 0; % Return 0 if no match is found
    end
end


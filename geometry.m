%% Spacecraft Geometry and Surface Properties

alpha_MLI = 0.1;
alpha_radiator = 0.15;
alpha_louver = 0.15;
epsilon_MLI = 0.1;
epsilon_radiator = 0.8;
epsilon_louver = 0.8;

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
    
    alphas(3, alphas(1, :)~=0) = alpha_louver;
    epsilons(3, epsilons(1, :)~=0) = epsilon_louver;
    
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

figure(1)
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

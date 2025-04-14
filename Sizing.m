clc
clear
%% Output file name to grab
name_output_file = 'VenusUranus_energies.mat';
if ispc
    addpath('Functions\', 'Input\', 'Output\');
    % Load Inc_Energies.m outputs (may need to edit filename)
    path_out_file = strcat('Output\',name_output_file);
elseif ismac || isunix
    addpath('Functions/', 'Input/', 'Output/');
    % Load Inc_Energies.m outputs (may need to edit filename)
    path_out_file = strcat('Output/',name_output_file);
else
    disp('Platform not supported!');
end
%% Geometry maker
% on the other hand can edit the geometry.m script and uncomment the
% following group comment but do not change the areas in the script unless
% you run the Inc_Energies.m script again

% runs geometry script
geometry
% resave geometry
save(path_out_file, "pat", "-append");

%% Setup
load(path_out_file)


%%%% NOTE: IF YOU WANT TO CHANGE THE BELOW STRUCT VALUES WITHOUT RUNNING
%%%% Inc_Energies.m, USE THE FUNCTION updatepatchdata FROM THE FUNCTIONS
%%%% FOLDER YOU CAN PASS THE UPDATED MATRIX OF ALPHAS AND EPSILONS IN AS
%%%% PARAMETERS. ON THE OTHER HAND, YOU CAN ALSO JUST SKIP THIS AND
%%%% WRITE OUT THE ALPHA AND EPSILON VALUES YOURSELF
% Grab struct data for ease of use
areas = pat.UserData.Areas;
alphaRad = pat.UserData.Alphas(1,:);
alphaMLI = pat.UserData.Alphas(2,:);
epsRad = pat.UserData.Epsilons(1,:);
epsMLI = pat.UserData.Epsilons(2,:);
alphaLouv = pat.UserData.Alphas(3, :);
epsLouv = pat.UserData.Alphas(3, :);


%% USER INPUTS
upperT = (80+273.15); % upper temperature bound [K]
lowerT = -20+273.15; % lower temperature bound [K]
uncertainty = 20; % uncertainty of temperature [K]
powE = 0.8*2000; % Thermal power input from internals [W] (estimated as 80% of RTG heat)
T_upper = upperT-uncertainty; % Target operating temp [K]
T_lower = lowerT+uncertainty; % Lower bound for non-operating case [K]

%% Run functions
% Hot Case radiator Sizing
%[areaRad,areaMLI,T_ss,err] = szRadiator(inc_S_avgs{1},inc_IR_avgs{1},inc_A_avgs{1},areas,alphaRad,alphaMLI,epsRad,epsMLI,powE*.90,T_upper,false);
%disp("Avg temp (hot case rad sizing): " + (T_ss-273.15) + " deg C")
% assuming 90% power operations required to maintain spacecraft (should
% change to actual value from power subsystem if they have it)

% Hard-coding
areaRad = [0 0 0 0 0 0 0 1.949051493563100 0 0];
areaMLI = areas-areaRad;

% cold case heater sizing using louver data
Q_htr = szHeater(inc_S_avgs{2},inc_IR_avgs{2},inc_A_avgs{2},areaRad,areas,alphaLouv,alphaMLI,epsLouv,epsMLI,powE,T_lower, T_upper);

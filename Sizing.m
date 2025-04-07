%% Setup
% Load Inc_Energies.m outputs (may need to edit filename)
load('VenusUranus_energies.mat')

% Grab struct data for ease of use
areas = pat.UserData.Areas;
alphaRad = pat.UserData.Alphas(1,:);
alphaMLI = pat.UserData.Alphas(2,:);
epsRad = pat.UserData.Epsilons(1,:);
epsMLI = pat.UserData.Epsilons(2,:);

if ispc
    addpath('Functions\', 'Input\', 'Output\');
elseif ismac || isunix
    addpath('Functions/', 'Input/', 'Output/');
else
    disp('Platform not supported!');
end

%% USER INPUTS
powE = 500; % Thermal power input from internals [W]
tgtTavg = (25+273.15)*ones(size(areas)); % Target operating temp for each face [K]
T_bound = 0+273.15; % Lower bound for non-operating case [K]

%% Run functions
for ii=1:length(inc_A_avgs)
    [areaRad{ii},areaMLI{ii},err{ii}] = szRadiator(inc_S_avgs{ii},inc_IR_avgs{ii},inc_A_avgs{ii},areas,alphaRad,alphaMLI,epsRad,epsMLI,powE,tgtTavg);
    Q_htr{ii} = szHeater(inc_S_avgs{ii},inc_IR_avgs{ii},inc_A_avgs{ii},areaRad{ii},areas,alphaRad,alphaMLI,epsRad,epsMLI,powE,T_bound);
end

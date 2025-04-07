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
    disp('Platform not supported!');
end
%% make the name of the file you saved
% planets you are working with

savename = "";
for planet = planets
    savename = strcat(savename, planet);
end
savename = strcat(savename, "_energies");

%% load the file with the energies, geometry, and thermal properties
load(savename)

%% Calculate the temperatures, size radiators, and heaters

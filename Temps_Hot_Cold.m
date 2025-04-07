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
% choose the planets you have saved data for
planets = ["Venus", "Uranus"];
% creates the path to the file
if ispc
    savename = "Output\";
elseif ismac || isunix
    savename = "Output/";
else
    error('Platform not supported!');
end
% creates the name of the file to load the data from along the path
for planet = planets
    savename = strcat(savename, planet);
end
savename = strcat(savename, "_energies");

%% load the file with the energies, geometry, and thermal properties
load(savename)

%% Calculate the temperatures, size radiators, and heaters

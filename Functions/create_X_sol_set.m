clear all
assume(symvar(assumptions), 'clear')
close all
clc

X_sol_set = table('RowNames',{'x_sol', 'y_sol', 'z_sol',...
    'conx', 'cony', 'conz'});

R_s = 1;
d = 2;

%{
% NOTE: special cases no longer apply because when the plane parallel to
% the surface and the plane of the spherical cap are the same, it is the
% same as saying there is no intersection because that is how it works out
in the decision matrix.

% special case 1
gamma = 0;
phi = 80*pi/180;
omega = 80*pi/180;
psi = acos(d.*cos(phi)./R_s);
[~, ~, ~, X_sol_set] =...
    X_intersect(d, R_s, gamma, phi, omega,...
        psi, X_sol_set);


% special case 2
gamma = pi;
phi = 80*pi/180;
omega = 80*pi/180;
psi = acos(d.*cos(phi)./R_s);
[~, ~, ~, X_sol_set] =...
    X_intersect(d, R_s, gamma, phi, omega,...
        psi, X_sol_set);


% special case 3
gamma = -pi;
phi = 80*pi/180;
omega = 80*pi/180;
psi = acos(d.*cos(phi)./R_s);
[~, ~, ~, X_sol_set] =...
    X_intersect(d, R_s, gamma, phi, omega,...
        psi, X_sol_set);

% special case 4
gamma = 0;
phi = 90*pi/180;
omega = 90*pi/180;
psi = acos(d.*cos(phi)./R_s);
[~, ~, ~, X_sol_set] =...
    X_intersect(d, R_s, gamma, phi, omega,...
        psi, X_sol_set);

% special case 5
gamma = 0;
phi = 100*pi/180;
omega = 100*pi/180;
psi = acos(d.*cos(phi)./R_s);
[~, ~, ~, X_sol_set] =...
    X_intersect(d, R_s, gamma, phi, omega,...
        psi, X_sol_set);

% special case 6
gamma = pi;
phi = 90*pi/180;
omega = 90*pi/180;
psi = acos(d.*cos(phi)./R_s);
[~, ~, ~, X_sol_set] =...
    X_intersect(d, R_s, gamma, phi, omega,...
        psi, X_sol_set);

% special case 7
gamma = -pi;
phi = 90*pi/180;
omega = 90*pi/180;
psi = acos(d.*cos(phi)./R_s);
[~, ~, ~, X_sol_set] =...
    X_intersect(d, R_s, gamma, phi, omega,...
        psi, X_sol_set);
%}

% set of values to go through different cases
gammas = -pi:45*pi/180:pi;
%gammas = gammas(gammas~=-pi & gammas~=pi & gammas~=0);
phis = 0:45*pi/180:pi;
omegas = 0:45*pi/180:pi;
omegas = omegas(omegas~=pi & omegas~=0);
psis = 0:45*pi/180:pi/2;
psis = psis(psis ~= 0);
for i = 1:length(phis)
    for j = 1:length(omegas)
        for k = 1:length(gammas)
            for l = 1:length(psis)
                [~, ~, ~, X_sol_set] =...
                    X_intersect(d, R_s, gammas(k), phis(i), omegas(j),...
                        psis(l), X_sol_set);
                %X_sol_set
            end
        end
    end
end
save("calc_sym_inter.mat", "X_sol_set")

%%
for i = 1:size(X_sol_set, 1)
    for j = 1:size(X_sol_set, 2)
        if length(X_sol_set{i, j}{1}) ~= 2
            warning("length not equal to 2 at %i, %i\n", i, j)
        end
    end
end


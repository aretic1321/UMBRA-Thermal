clc;
clear;
% eclipse time
% source: https://commons.erau.edu/cgi/viewcontent.cgi?article=1412&context=ijaaa
beta = 0;
h=500*10^3;
G = 6.6743*10^-11;
M = 5.9722*10^24;
R_e=6378.1366*10^3;
f_e = @(Beta) 1/180*acosd(sqrt(h^2+2*R_e*h)./((R_e+h)*cosd(Beta)));
Beta_star = asind(R_e/(R_e+h));
T=2*pi*sqrt((R_e+h)^3/(G*M));
if abs(beta) < Beta_star
    % fraction of time in sunlight
    f_s_eval = 1 - f_e(beta);
else
    % fraction of time in sunlight
    f_s_eval = 1 - 0;
end

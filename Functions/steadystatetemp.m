function T = steadystatetemp(areaRad,areaMLI,inc_S,inc_IR,inc_A,alphas,MLI_alpha,epsilons,MLI_epsilon,powE)
    % Stephan Boltzman Constant
    boltz = 5.67051*(10^-8);

    % Absorbed heat [W]
    absb_A = inc_A.*alphas.*areaRad + inc_A.*MLI_alpha.*areaMLI; % Albedo
    absb_S = inc_S.*alphas.*areaRad + inc_S.*MLI_alpha.*areaMLI; % Solar
    absb_IR = inc_IR.*epsilons.*areaRad + inc_IR.*MLI_epsilon.*areaMLI; % Earth IR
    
    % Total absorbed heat [W]
    Q_in = sum(powE) + sum(absb_IR + absb_A + absb_S); % Heat in [W]
    
    % Face temperature [K]
    T = ( Q_in./(sum(areaRad.*epsilons*boltz + areaMLI.*MLI_epsilon*boltz)) ).^(1/4);
end
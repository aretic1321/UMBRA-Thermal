    % FOR TESTING (COMMENT OUT)
    % clear
    % clc
    % close all
    % 
    % % Test values to match excel doc (Hot Case)
    % inc_S = [154.2 154.2 1333.4 0.0 154.2 154.2];
    % inc_IR = [0.0 224.6 69.9 70.0 69.7 69.7];
    % inc_A = [0.0 59.6 26.5 11.5 18.6 18.6];
    % areas = [1.00 1.00 1.00 1.00 1.00 1.00];
    % alphas = [0.15 0.15 0.15 0.15 0.15 0.15];
    % MLI_alpha = [0.10 0.10 0.10 0.10 0.10 0.10];
    % epsilons = [0.80 0.80 0.80 0.80 0.80 0.80];
    % MLI_epsilon = [0.02 0.02 0.02 0.02 0.02 0.02];
    % powE = [0 0 0 0 0 0];
    % tgtTavg = [68.32 -10.64 316.81 -55.39 85.07 85.07] + 273.15; % [K]
    % 
    % % Run test (rename 'szRadiator(...)' to 'szRadiatorTest(...)' )
    % [areaRad,areaMLI,err] = szRadiatorTest(inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg);

function [areaRad,areaMLI,err] = szRadiator(inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg)
    % RADIATOR SIZING FUNCTION
    
    % INPUTS
    % inc_S: Inclined solar energy/power (1 vector of values corresponding to the areas)
    % inc_IR: Inclined IR energy/power (1 vector of values corresponding to the areas)
    % inc_A: Inclined Albedo energy/power (1 vector of values corresponding to the areas)
    % areas: Total areas of each face (1 vector of values)
    % alphas: alpha value of each RADIATOR (1 vector of values)
    % MLI_alpha: alpha value of each MLI (1 vector of values)
    % epsilons: epsilon value of each RADIATOR (1 vector of values)
    % MLI_epsilon: epsilon value of each MLI (1 vector of values)
    % pow_E: applied electrical power
    % tgt_Tavg: Desired average temperature to maintain
    
    % OUTPUTS
    % areaRad: area of the radiator on each face (1 vector of values)
    % areaMLI: area of MLI on each face (1 vector of values)
    
    % Initialize with radiator area = 0 m2
    areaRad0 = 0*areas;
    
    % fsolve options
    optns = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','MaxIterations',1e6,'MaxFunctionEvaluations',1e6);

    % Iteratively solve for radiator areas [m2]
    [areaRad,err] = fsolve(@(x) objFunc(x,inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg),areaRad0,optns);
    
    % Results
    areaMLI = areas-areaRad; % MLI area result [m2]
end

% Objective function. Minimize T_err to zero
function T_err = objFunc(areaRad,inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg)
    % Stephan Boltzman Constant
    boltz = 5.67051*(10^-8); 
    
    % Update MLI areas
    areaMLI = areas-areaRad;
    
    % Absorbed heat [W]
    absb_A = inc_A.*alphas.*areaRad + inc_A.*MLI_alpha.*areaMLI; % Albedo
    absb_S = inc_S.*alphas.*areaRad + inc_S.*MLI_alpha.*areaMLI; % Solar
    absb_IR = inc_IR.*epsilons.*areaRad + inc_IR.*MLI_epsilon.*areaMLI; % Earth IR
    
    % Total absorbed heat [W]
    Q_in = powE + absb_IR + absb_A + absb_S; % Heat in [W]
    
    % Face temperature [K]
    T = ( Q_in./(areaRad.*epsilons*boltz + areaMLI.*MLI_epsilon*boltz) ).^(1/4);
    
    % Objective function: Difference between resulting temp and target temp [K]
    T_err = abs( (tgtTavg - T)./tgtTavg );
end
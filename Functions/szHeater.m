% FOR TESTING (COMMENT OUT)
    % clear
    % clc
    % close all
    % 
    % % Test values to match excel doc (Cold, non-operating Case)
    % areaRad = [0.00 0.74 0.00 0.00 0.00 0.00];
    % inc_S = [418.2 30.4 0.8 0.8 287.3 287.3];
    % inc_IR = [0.0 186.8 58.3 58.2 58.0 58.0];
    % inc_A = [0.0 79.1 24.5 24.5 24.6 24.6];
    % areas = [1.00 1.00 1.00 1.00 1.00 1.00];
    % alphas = [0.05 0.05 0.05 0.05 0.05 0.05];
    % MLI_alpha = [0.05 0.05 0.05 0.05 0.05 0.05];
    % epsilons = [0.80 0.80 0.80 0.80 0.80 0.80];
    % MLI_epsilon = [0.02 0.02 0.02 0.02 0.02 0.02];
    % powE = 0; % W
    % T_bound = 0+273.15; % K
    % 
    % % Run test (rename 'szHeater(...)' to 'szHeaterTest(...)' )
    % Q_htr = szHeaterTest(inc_S,inc_IR,inc_A,areaRad,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,T_bound);

function Q_htr = szHeater(inc_S,inc_IR,inc_A,areaRad,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,T_bound)
    % HEATER SIZING FUNCTION
    
    % INPUTS
    % inc_S: Inclined solar energy/power (1 vector of values corresponding to the areas)
    % inc_IR: Inclined IR energy/power (1 vector of values corresponding to the areas)
    % inc_A: Inclined Albedo energy/power (1 vector of values corresponding to the areas)
    % areas: Total areas of each face (1 vector of values)
    % areaRad: Areas of each RADIATOR (1 vector of values)
    % alphas: alpha value of each RADIATOR (1 vector of values)
    % MLI_alpha: alpha value of each MLI (1 vector of values)
    % epsilons: epsilon value of each RADIATOR (1 vector of values)
    % MLI_epsilon: epsilon value of each MLI (1 vector of values)
    % T_bound: Desired average temperature to maintain (1 value) [K]
    
    % OUTPUTS
    % Q_htr: Heat contribution from heater req'd to meet target temps [W]
    
    % Average powE in to each side
    powE = powE*ones(size(areas))/length(areas);

    % Stephan Boltzman Constant
    boltz = 5.67051*(10^-8); 
    
    % Update MLI areas
    areaMLI = areas-areaRad;
    
    % Absorbed heat [W] *empty MLI alpha and epsilon variables
    absb_A = inc_A.*alphas.*areaRad + inc_A.*MLI_alpha.*areaMLI; % Albedo
    absb_S = inc_S.*alphas.*areaRad + inc_S.*MLI_alpha.*areaMLI; % Solar
    absb_IR = inc_IR.*epsilons.*areaRad + inc_IR.*MLI_epsilon.*areaMLI; % Earth IR
    
    % Total absorbed heat [W] (per face)
    Q_in = powE + absb_IR + absb_A + absb_S;
    
    % Face temperature [K] (per face)
    %T = ( Q_in./(areaRad.*epsilons*boltz + areaMLI.*MLI_epsilon*boltz) ).^(1/4);

    % Total heat lost to surroundings [W] (sum)
    Q_out = sum(boltz*(T_bound.^4).*(areaRad.*epsilons + areaMLI.*MLI_epsilon));
    
    % Heater Contribution Required [W]
    Q_htr = Q_out-sum(Q_in);
end
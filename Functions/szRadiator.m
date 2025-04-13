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

    function [areaR,areaMLI,T,err] = szRadiator(inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg, minarea)
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
    % minarea: truth value of whether to find the minimu area or maximum
    
    % OUTPUTS
    % areaRad: area of the radiator on each face (1 vector of values)
    % areaMLI: area of MLI on each face (1 vector of values)
    
    % Initialize with radiator area = 0 m2
    areaRad0 = zeros([1 length(areas(alphas~=0))]);
    areaRad0 = 0*areas(alphas~=0);
    
    %{
    % fsolve options
    optns = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt','MaxIterations',1e6,'MaxFunctionEvaluations',1e6, 'Display', 'iter');

    % Iteratively solve for radiator areas [m2]
    [areaRad,err] = fsolve(@(x) objFunc(x,inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg),areaRad0,optns);
    %}

    % fmincon options
    optns = optimoptions(@fmincon,'Algorithm', 'sqp','MaxIterations',1e6,'MaxFunctionEvaluations',1e6);
    f = @(x) minfun(x,inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg,minarea);
    problem = createOptimProblem('fmincon','x0',areaRad0,...
        'objective',f,'lb',zeros([length(areas(alphas~=0)) 1]),'ub',areas(alphas~=0), 'options', optns);
    [areaRad, fout] = run(GlobalSearch,problem);
    %{
    [areaRad, fout] = fmincon(,...
        areaRad0',... -1*diag(ones([length(areas(alphas~=0)) 1])), zeros([length(areas(alphas~=0)) 1]));
        [], [], [], [], zeros([length(areas(alphas~=0)) 1]), areas(alphas~=0));
    %}
    err = abs(fout) - sum(areaRad);
    % Results
    areaR = zeros([1 length(areas)]);
    areaR(alphas~=0) = areaRad';

    areaMLI = areas-areaR; % MLI area result [m2]
    
    T = steadystatetemp(areaR, areaMLI,inc_S,inc_IR,inc_A,alphas,MLI_alpha,epsilons,MLI_epsilon,powE);
    if err > 1e-4
        error("Target temperature is infeasible for the current " +...
            "MLI and Radiator Properties. Achievable temperature is " +...
            "%d K, while target is %d K.", T, tgtTavg)
    end
end

%{
% Objective function. Minimize T_err to zero
function T_err = objFunc(areaRad,inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg)
    
    areaR = zeros([1 length(areas)]);
    areaR(alphas~=0) = areaRad;

    % Stephan Boltzman Constant
    boltz = 5.67051*(10^-8); 
    
    % Update MLI areas
    areaMLI = areas-areaR;
    
    % Absorbed heat [W]
    absb_A = inc_A.*alphas.*areaR + inc_A.*MLI_alpha.*areaMLI; % Albedo
    absb_S = inc_S.*alphas.*areaR + inc_S.*MLI_alpha.*areaMLI; % Solar
    absb_IR = inc_IR.*epsilons.*areaR + inc_IR.*MLI_epsilon.*areaMLI; % Earth IR
    
    % Total absorbed heat [W]
    Q_in = powE + sum(absb_IR + absb_A + absb_S); % Heat in [W]
    
    % Face temperature [K]
    T = ( Q_in./(sum(areaR.*epsilons*boltz + areaMLI.*MLI_epsilon*boltz)) ).^(1/4);
    
    % Objective function: Difference between resulting temp and target temp [K]
    T_err = abs( (tgtTavg - T)./tgtTavg );
end
%}

function out = minfun(areaRad,inc_S,inc_IR,inc_A,areas,alphas,MLI_alpha,epsilons,MLI_epsilon,powE,tgtTavg, minarea)
    areaR = zeros([1 length(areas)]);
    areaR(alphas~=0) = areaRad';
    
    % Update MLI areas
    areaMLI = areas-areaR;
    out = abs(steadystatetemp(areaR, areaMLI,inc_S,inc_IR,inc_A,alphas,MLI_alpha,epsilons,MLI_epsilon,powE)...
            - tgtTavg);
    if minarea
        out = out+sum(areaRad);
    else
        out = out-sum(areaRad);
    end
end


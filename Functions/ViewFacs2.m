function [Fs, varargout] = ViewFacs2(R_s, ds, phis, omegas, gammas, varargin)
    % Finds the view factors
    % using what was found in 
    % https://www.sciencedirect.com/science/article/pii/S001793102033413X?fr=RR-2&ref=pdf_download&rr=8f4ba51f696f9123#fig0003
    % Hoewever, there were some errors in the article and have been
    % fixed in this program.
    % properly adjusted.
    % required inputs:
    %   R_s: Radius of full sphere that the cap is part of 
    %   ds: distance of the infinitesimal sphere to the center of the
    %   sphere the cap is a part of
    %   phis: central angles between center of spherical cap to -z direction  
    %   omegas: angles between the z axis and the normal vector
    %   gammas: angles between the +x direction the normal surface
    % parametric inputs:
    %   normals: infinitesimal surface normal
    %   psi: central angle between center of spherical cap to its edge
    %   X_sol_sets: table of intersection solutions for the plane of the
    %   plate, projection plane, and sphere
    % outputs:
    %   Fs: view factors
    %   varargout: table of intersection solutions for the plane of the
    %   plate, projection plane, and sphere
    % note: psi should specified if light can't be assumed to only cover
    % half of a planet
    % note: infinitesimal surface is on the z-axis in the negative
    % direction relative to the sphere

    p = inputParser;
    addParameter(p, 'psi', 90);
    addParameter(p, 'debug', false);
    addParameter(p, 'X_sol_set',...
        table('RowNames',{'x_sol', 'y_sol', 'z_sol', 'conx', 'cony',...
        'conz'}), @(T) isempty(T) || istable(T));

    parse(p, varargin{:});

    debug = p.Results.debug;
    X_sol_set = p.Results.X_sol_set;

    Fs = zeros(1, length(omegas));

    if isempty(X_sol_set) && exist("calc_sym_inter.mat", 'file')
        load("calc_sym_inter.mat","X_sol_set")
    else
        X_sol_set = table('RowNames',{'x_sol', 'y_sol', 'z_sol', 'conx', 'cony',...
            'conz'});
    end
    for i = 1:length(omegas)
        [Fs(i), X_sol_set] = ViewFac(...
            R_s, ds(i), phis(i), omegas(i), gammas(i),...
            'debug', debug, 'X_sol_set', X_sol_set);
    end
    
    varargout{1} = X_sol_set;
end
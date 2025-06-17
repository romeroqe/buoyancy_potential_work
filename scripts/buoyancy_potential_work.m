% Function to calculate the WORK DONE BY THE BUOYANCY FORCE (WB) in
% vertically displacing a water parcel from different depths to a reference
% depth zref.
function [WB, z] = buoyancy_potential_work(rho, z, zref)
    % WB quantifies the water column's vertical homogeneity in terms of the work
    % done by the buoyancy force in vertically displacing a water parcel under 
    % static instability conditions.
    %
    % Inputs:
    % * rho: Sigma-0. Potential Density Anomaly profile referred to 0 dbar (kg m^-3)
    % * z: Depths in negative values (m)
    % * zref: Reference depth (default -10 m). WB(zref)=0 by definition.
    % The data must be a column vector and be arranged from the greatest to the shallowest depth.
    %
    % Outputs:
    % * WB: Vertical profile of work done by buoyancy (J m^-3)
    % * z: Depths (m)

    % If zref is not provided, it is assigned to its default value
    if nargin < 3
        zref = -10;
    end

	% If there is no data at zref, it is interpolated
    if ~any(z == zref)
        % Interpolate rho at zref
        rho_zref = interp1(z, rho, zref, 'linear', 'extrap');

        % Construct the vectors z and rho, adding zref and rho_zref
        z_extended = [zref; z];
        rho_extended = [rho_zref; rho];

        % Sort vectors z and rho
        [z,order] = sort(z_extended, 'ascend');
        rho = rho_extended(order);
    end

    % Defining variables
    g = 9.81; % Acceleration of gravity (m s^-2)
    nz = length(z); % Number of data in the vertical
    WB = NaN(nz,1); % Work done by buoyancy (J m^-3)
    izref = find(z==zref); % Index of zref in the vector z
    dz = diff(z); % Differences between adjacent elements of z (m), with positive units

    % Calculation of the vertical profile of WB
    for it = 1:nz % Loop for every depth
        zint = z(it);
        if isnan(zint) % If the depth of analysis is nan, WB is also nan
            WB(it) = NaN;
        else
            rho_int = rho(it);

            % Computation of WB
            S = NaN(nz,1);
            S(it) = 0;

            for i = 1:1:it-1
                sr = rho(i:it-1) + rho(i+1:it);
                sdz = dz(i:it-1);
                S(i) = -0.5 * sum(sr .* sdz);
            end
            for i = it+1:1:nz
                sr = rho(it:i-1) + rho(it+1:i);
                sdz = dz(it:i-1);
                S(i) = 0.5 * sum(sr .* sdz);
            end

            V = -g * (z - zint) * rho_int + g * S; % Buoyancy potential (J m^-3)
            V = -V;  % Change the sign of the buoyancy potential
            WB(it) = V(izref)-V(it); % WB at the depth of analysis (J m^-3)
        end
    end
end
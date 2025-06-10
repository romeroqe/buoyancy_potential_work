function [WB, z] = buoyancy_potential_work(rho, z, zint)
    % Calculates the WORK DONE BY BUOYANCY (WB)
    % WB quantifies the water column's vertical homogeneity in terms of the work
    % done by the buoyancy force in vertically displacing a water parcel under 
    % static instability conditions.
    %
    % Inputs:
    % - rho: Sigma-0 Potential Density Anomaly profile (kg/m^3), column vector
    % - z: Depths (m), negative values, deepest to shallowest (also column vector)
    % - zint: Reference depth (default -10 m)
    %
    % Outputs:
    % - WB: Work done by buoyancy profile (J/m^3)
    % - z: Depths (m)

    if nargin < 3
        zint = -10;
    end

	% A shallower depth than that of interest is required
    if z(1) < zint
        WB = [];
        z = [];
        return;
    end

	% If there is no data at zint, it is interpolated
    if ~any(z == zint)
        rho_interp = interp1(z, rho, zint, 'linear');
        i10 = find(z <= zint, 1);
        z = [z(1:i10-1); zint; z(i10:end)];
        rho = [rho(1:i10-1); rho_interp; rho(i10:end)];
    end

    nz = length(z); % Amount of data in the vertical
    g = 9.81; % Acceleration of gravity
    WB = NaN(nz,1);

    izint = find(z == zint, 1);
    dz = diff(z); % Depth vector spacing (m) with positive units

    % Calculate buoyancy work for various depths of interest
    for it = 1:nz
        zint = z(it);
        if isnan(zint)
            % If the depth of interest is nan, the work is also nan
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

            V = -g * (z - zint) * rho_int + g * S; % Buoyancy potential (J*m^-3)
            V = -V;  % Change the sign of the potential to make it more intuitive
            WB(it) = V(izint); % Work done by buoyancy (J*m^-3)
        end
    end
end

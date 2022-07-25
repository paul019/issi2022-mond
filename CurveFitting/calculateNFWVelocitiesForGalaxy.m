function [r, NFWVelocities] = calculateNFWVelocitiesForGalaxy(rotationCurveData,p0,R_S)

% Define constants:
G = 6.6743e-2;                      % in km^3 * kg^(-1) * s^(-2)

r = rotationCurveData(:,9);         % in km
Vbaryon = rotationCurveData(:,11);  % in km/s

V_DM = sqrt( G * r.^(-1) * 4 * pi * p0 * R_S^3 .* ( log( (R_S + r) ./ R_S ) + R_S ./ (R_S + r) - 1 ) );

NFWVelocities = sqrt( Vbaryon.^2 + V_DM.^2 );

end


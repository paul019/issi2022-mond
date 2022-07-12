function [r, MONDVelocities] = calculateMONDVelocitiesForGalaxy(rotationCurveData,a0,interpolationFunction)

% Default value for interpolationFunction:
if nargin < 3
    interpolationFunction = 'linear';
end

r = rotationCurveData(:,9);         % in km
Vbaryon = rotationCurveData(:,11);  % in km/s

% Calculate MOND velocity in km/s:
switch interpolationFunction
    case 'linear'
        % Linear approach:
        MONDVelocities = nthroot(abs((Vbaryon.^2).*r*a0),4);

    case 'rar'
        % RAR approach:
        MONDVelocities=sqrt((real(Vbaryon).^2)./(1-exp(-sqrt((real(Vbaryon).^2)./(r * a0)))));

    % PUT IN MORE INTERPOLATION FUNCTIONS

    otherwise
        % Linear approach:
        MONDVelocities = nthroot(abs((Vbaryon.^2).*r*a0),4);
end

end


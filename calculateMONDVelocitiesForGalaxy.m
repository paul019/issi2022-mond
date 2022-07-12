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

    case 'simple'
        % Simple interpolation curve approach:
        MONDVelocities=sqrt(((Vbaryon.^2) + sqrt((real(Vbaryon).^4) + 4*(real(Vbaryon).^2).* a0.*r))/2);

    case 'standard'
        % Standard interpolation curve approach:
        MONDVelocities=nthroot(((Vbaryon.^4) + sqrt((Vbaryon.^8) + 4*(Vbaryon.^4).* (a0^2).*(r).^2))/2 ,4);

        % PUT IN MORE INTERPOLATION FUNCTIONS

    otherwise
        % Linear approach:
        MONDVelocities = nthroot(abs((Vbaryon.^2).*r*a0),4);
end

end


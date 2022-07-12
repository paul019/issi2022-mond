function [r, MONDVelocities] = calculateMONDVelocitiesForGalaxy(rotationCurveData,a0)

r = rotationCurveData(:,9);         % in km
Vbaryon = rotationCurveData(:,11);  % in km/s

% Calculate MOND velocity in km/s:

    % Linear approach:
    %MONDVelocities = nthroot(abs((Vbaryon.^2).*r*a0),4);
    
    % RAR approach:
    MONDVelocities=sqrt((real(Vbaryon).^2)./(1-exp(-sqrt((real(Vbaryon).^2)./(r * a0)))));

end


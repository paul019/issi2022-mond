function [r, MONDVelocities] = calculateMONDVelocitiesForGalaxy(rotationCurveData,a0)

r = rotationCurveData(:,9);
Vbaryon = rotationCurveData(:,11);

% MOND velocity in km/s

%MONDVelocities = nthroot(abs((Vbaryon.^2).*r*a0),4);
MONDVelocities=sqrt((real(Vbaryon).^2)./(1-exp(-sqrt((real(Vbaryon).^2)./(r * a0)))));

end


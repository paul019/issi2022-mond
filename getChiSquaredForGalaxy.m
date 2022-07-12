function chiSquared = getChiSquaredForGalaxy(rotationCurveData,a0)

% Calculate MOND velocities:
[~, MONDVelocities] = calculateMONDVelocitiesForGalaxy(rotationCurveData,a0);

% Get observed velocities:
observedVelocities = rotationCurveData(:,2);        % in km/s

% Get the errors on the observed velocities:
observedVelocitiesError = rotationCurveData(:,3);   % in km/s

% Calculate the mean squared weighted deviation:
chiSquared = sum( ( (observedVelocities - MONDVelocities) ./ observedVelocitiesError ).^2 );
    % source: https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic

end


function chiSquared = getChiSquaredForGalaxy(dataOrName,a0)

if ischar(dataOrName)
    data = prepareGalaxyForMOND(dataOrName,0.5,0.5);
else
    data = dataOrName;
end

[~, MONDVelocities] = calculateMONDVelocitiesForGalaxy(data,a0);
observedVelocities = data(:,2);
observedVelocitiesError = data(:,3);

chiSquared = sum(((observedVelocities - MONDVelocities).^2 ./ observedVelocitiesError).^2);

end


function [a0Values, chiSquaredVector] = findBestA0ForGalaxy(name,a0Min,a0Step,a0Max)

galaxyData = prepareGalaxyForMOND(name,0.5,0.7);

a0Values = transpose(a0Min:a0Step:a0Max);
chiSquaredVector = zeros(length(a0Values), 1);

for ii = 1:length(a0Values)
    progressInPercent = (ii-1) / length(a0Values) * 100;
    fprintf('%f %%\n', progressInPercent);
    chiSquaredVector(ii) = getChiSquaredForGalaxy(galaxyData,a0Values(ii));
end

progressInPercent = 100;
fprintf('%f %%\n', progressInPercent);

figure
plot(a0Values, chiSquaredVector);

end


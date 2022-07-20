function plotChiSquaredHistogram(galaxyNames, galaxyFittingData, worstChiSquaredReduced)

% Sort galaxies by chiSquaredReduced:
[galaxyNames, galaxyFittingData] = sortGalaxiesByChiSquared(galaxyNames, galaxyFittingData, false);

% Iterate through the sorted galaxies and get rid of bad ones:
jj = 1;
while galaxyFittingData{jj}.chiSquaredReduced_general <= worstChiSquaredReduced
    jj = jj + 1;
end
galaxyNames = galaxyNames{1:jj};

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

% Create a vector with the minimum chi squared for every galaxy:
chiSquaredReduced = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    chiSquaredReduced(jj) = galaxyFittingData{jj}.chiSquaredReduced_general;
end

figure

histogram(chiSquaredReduced, 200);

set(gca,'FontSize',15);
xlabel '\chi_\nu^2';
xlim([0 10]);
ylabel 'Number of galaxies';
title('Histogram of all galaxies by \chi_\nu^2');
grid on;

end


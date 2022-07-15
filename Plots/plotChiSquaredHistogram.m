function plotChiSquaredHistogram(galaxyNames, galaxyFittingData)

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

% Create a vector with the minimum chi squared for every galaxy:
chiSquaredReducedMin = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    chiSquaredReducedMin(jj) = galaxyFittingData{jj}.chiSquaredReducedMin;
end

figure

histogram(chiSquaredReducedMin, 200);

set(gca,'FontSize',15);
xlabel '\chi_\nu^2';
ylabel 'Number of galaxies';
title('Histogram of all galaxies by \chi_\nu^2');
grid on;

end


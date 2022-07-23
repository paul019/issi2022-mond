function plotChiSquaredHistogram(galaxyNames, galaxyFittingData, worstChiSquaredReduced, standaloneflag, title_, subtitle1)

if nargin < 5
    title_ = 'Histogram of all galaxies by \chi_\nu^2';
end

if nargin < 6
    subtitle1 = '';
end

totalNumOfGalaxies = length(galaxyNames);

% Sort galaxies by chiSquaredReduced:
[~, galaxyFittingData] = sortGalaxiesByChiSquared(galaxyNames, galaxyFittingData, false);

% Iterate through the sorted galaxies and get rid of bad ones:
numOfGalaxies = 0;
while galaxyFittingData{numOfGalaxies+1}.chiSquaredReduced_general <= worstChiSquaredReduced
    numOfGalaxies = numOfGalaxies + 1;

    if numOfGalaxies == totalNumOfGalaxies
        break
    end
end

% Create a vector with the minimum chi squared for every galaxy:
chiSquaredReduced = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    chiSquaredReduced(jj) = galaxyFittingData{jj}.chiSquaredReduced_general;
end

% Create second subtitle:
subtitle2 = sprintf('(%d of %d galaxies shown)', numOfGalaxies, totalNumOfGalaxies);

if standaloneflag
    figure
end

histogram(chiSquaredReduced, worstChiSquaredReduced);

set(gca,'FontSize',15);
xlabel '\chi_\nu^2';
%xlim([0 worstChiSquaredReduced]);
%ylim([0 20]);
ylabel 'Number of galaxies';
title(title_);
subtitle({subtitle1,subtitle2});
grid on;

end


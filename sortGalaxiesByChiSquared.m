function [galaxyNames, galaxyFittingData] = sortGalaxiesByChiSquared(galaxyNames, galaxyFittingData, printflag)

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

% Cut off the last entry representing the average of all galaxies:
if length(galaxyFittingData) > numOfGalaxies
    galaxyFittingData = galaxyFittingData(1:end-1);
end

galaxyFittingData_table = struct2table(cell2mat(galaxyFittingData));
[~, order] = sortrows(galaxyFittingData_table, 'chiSquaredReduced_general');

galaxyNames = galaxyNames(order);
galaxyFittingData = galaxyFittingData(order);

if printflag
    for ii = 1:numOfGalaxies
        fprintf('#%d %s: chiSquaredReducedMin = %d\n', ii, galaxyNames{ii}, galaxyFittingData{ii}.chiSquaredReduced_forBestOverallA0);
    end
end

end


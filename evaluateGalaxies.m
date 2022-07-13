function [galaxyNames, galaxyFittingData] = evaluateGalaxies(a0Min,a0Step,a0Max,interpolationFunction,galaxyNames)

% Default value for interpolationFunction:
if nargin < 4
    interpolationFunction = 'linear';
end

% Import the names and metadata of all galaxies from the SPARC database:
[allGalaxyNames,~] = ReadLelliC;

% Deal with the input argument "galaxyNames":
if nargin < 5
    galaxyNames = allGalaxyNames;
else
    myGalaxyNames = galaxyNames;

    galaxyNames = {};
    index = 1;

    for ii = 1:length(allGalaxyNames)
        for jj = 1:length(myGalaxyNames)
            if strcmp(allGalaxyNames{ii}, myGalaxyNames{jj})
                galaxyNames{index} = allGalaxyNames{ii};
                index = index + 1;
                break
            end
        end
    end
end

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

% Create a cell array for all galaxies:
galaxyFittingData = cell(numOfGalaxies + 1, 1);

% Iterate through all galaxies in order to get rotation curve data:
for ii = 1:numOfGalaxies
    galaxyFittingData{ii}.rotationCurveData = prepareGalaxyForMOND(galaxyNames{ii},0.5,0.7);
end

% Create a vector with all values of a0 we want to try out:
a0Values = transpose(a0Min:a0Step:a0Max);

% Create a matrix with one cell for every combination of a value of a0 and
% a galaxy:
chiSquared = zeros(length(a0Values), numOfGalaxies + 1);

% Iterate over all galaxies:
for jj = 1:numOfGalaxies

    % Iterate over all values of a0:
    for ii = 1:length(a0Values)

        % Calculating chi squared for the current combination of a0 and a
        % galaxy:
        chiSquared(ii,jj) = getChiSquaredForGalaxy(galaxyFittingData{jj}.rotationCurveData, a0Values(ii), interpolationFunction);
    end

    galaxyFittingData{jj}.a0Values = a0Values;

    % Each galaxy gets a vector with chi squared for every value of a0:
    galaxyFittingData{jj}.chiSquared = chiSquared(:,jj);

    % Get the minimum chi squared for this galaxy:
    galaxyFittingData{jj}.chiSquaredMin = min(chiSquared(:,jj));

    % Find the index of the minimum chi squared in the vector
    % "chiSquared":
    minPos = find(galaxyFittingData{jj}.chiSquared == galaxyFittingData{jj}.chiSquaredMin);
    galaxyFittingData{jj}.bestA0_index = minPos(1);

    % Find the value of a0, that corresponds to the minimum value of chi
    % squared; this value of a0 is the best a0 value:
    galaxyFittingData{jj}.bestA0 = a0Values(minPos(1));

    % Print galaxy data to console:
    %fprintf('%s: a_0 = %d km/s^2; chi^2 = %d km^2/s^2\n', galaxies{jj}.name, galaxies{jj}.bestA0, galaxies{jj}.chiSquaredMin);
end

% The last entry in the cell array "galaxies" represents the average of all
% galaxies:
jj = numOfGalaxies + 1;

galaxyFittingData{jj}.a0Values = a0Values;

% Calculate the mean chi squared for every value of a0 by averaging over
% all galaxies:
galaxyFittingData{jj}.chiSquared = mean(chiSquared,2);

% Get the minimum mean chi squared:
galaxyFittingData{jj}.chiSquaredMin = min(galaxyFittingData{jj}.chiSquared);

% Find the index of the minimum mean chi squared in the vector
% "chiSquaredMean":
minPos = find(galaxyFittingData{jj}.chiSquared == galaxyFittingData{jj}.chiSquaredMin);
galaxyFittingData{jj}.bestA0_index = minPos(1);

% Find the best value of a0 for all galaxies:
galaxyFittingData{jj}.bestA0 = a0Values(minPos(1));

% Print general findings to the console:
%fprintf('a_0 = %d; chi^2 = %d km^2/s^2\n', galaxies{jj}.bestA0, galaxies{jj}.chiSquaredMin);

end


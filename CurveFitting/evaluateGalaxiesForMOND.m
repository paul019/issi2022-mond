function [galaxyNames, galaxyFittingData] = evaluateGalaxiesForMOND(a0Min,a0Step,a0Max,interpolationFunctionId,printflag,galaxyNames,a0)

% Define constants:
MtoLdisk = 0.5;
MtoLbulge = 0.7;

% Default value for interpolationFunction:
if nargin < 5
    interpolationFunctionId = 'linear';
end

% Import the names and metadata of all galaxies from the SPARC database:
[allGalaxyNames,~] = ReadLelliC;

% Deal with the input argument "galaxyNames":
if nargin < 6
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

% Keep track of the total number of datapoints:
totalNumberOfDatapoints = 0;

%--------------------------------------------------------------------------

% Iterate through all galaxies in order to get rotation curve data:
for ii = 1:numOfGalaxies
    galaxyFittingData{ii}.rotationCurveData = prepareGalaxyRotationCurveData(galaxyNames{ii},MtoLdisk,MtoLbulge);
    totalNumberOfDatapoints = totalNumberOfDatapoints + length(galaxyFittingData{ii}.rotationCurveData);

    galaxyFittingData{ii}.typeOfFit = 'MOND';
end

% Print status update:
if printflag
    fprintf('\nEvaluating %d datapoints from %d galaxies.\n', totalNumberOfDatapoints, numOfGalaxies);
end

% Create a vector with all values of a0 we want to try out:
a0Values = transpose(a0Min:a0Step:a0Max);

% Create a matrix with one cell for every combination of a value of a0 and
% a galaxy:
chiSquared = zeros(length(a0Values), numOfGalaxies + 1);

%--------------------------------------------------------------------------

% Iterate over all galaxies:
for jj = 1:numOfGalaxies

    % Iterate over all values of a0:
    for ii = 1:length(a0Values)

        % Calculating chi squared for the current combination of a0 and a
        % galaxy:
        chiSquared(ii,jj) = getChiSquaredForGalaxy(galaxyFittingData{jj}.rotationCurveData, a0Values(ii), interpolationFunctionId);
    end

    galaxyFittingData{jj}.a0Values = a0Values;
    galaxyFittingData{jj}.intFctId = interpolationFunctionId;

    % Calculate the degrees of freedom:
    galaxyFittingData{jj}.numberOfDatapoints = length(galaxyFittingData{jj}.rotationCurveData);
    galaxyFittingData{jj}.degreesOfFreedom = length(galaxyFittingData{jj}.rotationCurveData) - 2;

    % Each galaxy gets a vector with chi squared for every value of a0:
    galaxyFittingData{jj}.chiSquared = chiSquared(:,jj);
    galaxyFittingData{jj}.chiSquaredReduced = chiSquared(:,jj) / galaxyFittingData{jj}.degreesOfFreedom;

    % Get the minimum chi squared for this galaxy:
    galaxyFittingData{jj}.chiSquaredMin = min(chiSquared(:,jj));
    galaxyFittingData{jj}.chiSquaredReducedMin = min(chiSquared(:,jj) / galaxyFittingData{jj}.degreesOfFreedom);

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

%--------------------------------------------------------------------------

% The last entry in the cell array "galaxies" represents the average of all
% galaxies:
jj = numOfGalaxies + 1;

galaxyFittingData{jj}.a0Values = a0Values;
galaxyFittingData{jj}.intFctId = interpolationFunctionId;

% Calculate the degrees of freedom:
galaxyFittingData{jj}.numberOfDatapoints = totalNumberOfDatapoints;
galaxyFittingData{jj}.degreesOfFreedom = totalNumberOfDatapoints - 2;

% Calculate the mean chi squared for every value of a0 by averaging over
% all galaxies:
galaxyFittingData{jj}.chiSquared = sum(chiSquared,2);
galaxyFittingData{jj}.chiSquaredReduced = sum(chiSquared,2) / galaxyFittingData{jj}.degreesOfFreedom;

% Get the minimum mean chi squared:
galaxyFittingData{jj}.chiSquaredMin = min(galaxyFittingData{jj}.chiSquared);
galaxyFittingData{jj}.chiSquaredReducedMin = min(galaxyFittingData{jj}.chiSquaredReduced);

galaxyFittingData{jj}.chiSquared_forBestOverallA0 = galaxyFittingData{jj}.chiSquaredMin;
galaxyFittingData{jj}.chiSquaredReduced_forBestOverallA0 = galaxyFittingData{jj}.chiSquaredReducedMin;

galaxyFittingData{jj}.chiSquared_general = galaxyFittingData{jj}.chiSquaredMin;
galaxyFittingData{jj}.chiSquaredReduced_general = galaxyFittingData{jj}.chiSquaredReducedMin;

% Find the index of the minimum mean chi squared in the vector
% "chiSquaredMean":
minPos = find(galaxyFittingData{jj}.chiSquared == galaxyFittingData{jj}.chiSquaredMin);
galaxyFittingData{jj}.bestA0_index = minPos(1);

if nargin < 7
    % Find the best value of a0 for all galaxies:
    galaxyFittingData{jj}.bestA0 = a0Values(minPos(1));
else
    % Use the value of a0 passed into the 
    galaxyFittingData{jj}.bestA0 = a0;
    galaxyFittingData{jj}.bestA0_index = find(a0Values == a0);
end

% Print general findings to the console:
if printflag
    fprintf('All galaxies: a_0 = %d; chi_v^2 = %d\n\n', galaxyFittingData{jj}.bestA0, galaxyFittingData{jj}.chiSquaredReducedMin);
end

bestA0Overall = galaxyFittingData{end}.bestA0;
bestA0Overall_index = galaxyFittingData{end}.bestA0_index;

%--------------------------------------------------------------------------

% Iterate over all galaxies:
for jj = 1:numOfGalaxies
    galaxyFittingData{jj}.chiSquared_forBestOverallA0 = galaxyFittingData{jj}.chiSquared(bestA0Overall_index);
    galaxyFittingData{jj}.chiSquaredReduced_forBestOverallA0 = galaxyFittingData{jj}.chiSquaredReduced(bestA0Overall_index);

    galaxyFittingData{jj}.chiSquared_general = galaxyFittingData{jj}.chiSquared_forBestOverallA0;
    galaxyFittingData{jj}.chiSquaredReduced_general = galaxyFittingData{jj}.chiSquaredReduced_forBestOverallA0;

    % Save galaxy data as human-readable string:
    galaxyFittingData{jj}.dataString = sprintf('%s: chi_v^2 = %d', galaxyNames{jj}, galaxyFittingData{jj}.chiSquaredReduced_general);

    % Save galaxy data as mond fit (input to
    % "plotGalaxyRotationCurveWithFits"):
    galaxyFittingData{jj}.mondFit = struct('intFctId',interpolationFunctionId,'a0',bestA0Overall,'chiSquaredReduced',galaxyFittingData{jj}.chiSquaredReduced_general);

    % Print galaxy data to console:
    if printflag
        fprintf('%s\n', galaxyFittingData{jj}.dataString);
    end
end

%--------------------------------------------------------------------------

bestChiSquaredReduced = 10^10;
bestGalaxyName = '';

for jj = 1:numOfGalaxies
    galaxy_chiSquaredReduced = galaxyFittingData{jj}.chiSquaredReduced_forBestOverallA0;

    if galaxy_chiSquaredReduced < bestChiSquaredReduced
        bestChiSquaredReduced = galaxy_chiSquaredReduced;
        bestGalaxyName = galaxyNames{jj};
    end
end

if printflag
    fprintf('\nBest galaxy for %s: %s (a_0 = %s m/s^2, chi_v^2 = %s)\n\n', interpolationFunctionId, bestGalaxyName, num2str(bestA0Overall*1e3), num2str(bestChiSquaredReduced));
end

end


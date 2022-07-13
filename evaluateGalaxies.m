function [a0Values, chiSquared, galaxies] = evaluateGalaxies(a0Min,a0Step,a0Max,interpolationFunction,galaxyNames)

% Default value for interpolationFunction:
if nargin < 4
    interpolationFunction = 'linear';
end

% Import the names and metadata of all galaxies from the SPARC database:
[allGalaxyNames,allGalaxyData] = ReadLelliC;

% Deal with the input argument "galaxyNames":
if nargin < 5
    galaxyNames = allGalaxyNames;
    galaxyData = allGalaxyData;
else
    myGalaxyNames = galaxyNames;

    galaxyNames = {};
    galaxyData = {};
    index = 1;

    for ii = 1:length(allGalaxyNames)
        for jj = 1:length(myGalaxyNames)
            if strcmp(allGalaxyNames{ii}, myGalaxyNames{jj})
                galaxyNames{index} = allGalaxyNames{ii};
                galaxyData{index} = allGalaxyData{ii};
                index = index + 1;
                break
            end
        end
    end
end

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

% Create a cell array for all galaxies:
galaxies = cell(numOfGalaxies,1);

% Iterate through all galaxies in order to get galaxy name, hubble type and
% rotation curve data:
for ii = 1:numOfGalaxies
    galaxies{ii}.name = galaxyNames{ii};
    galaxies{ii}.hubbleType = galaxyData{ii}(1);
    galaxies{ii}.rotationCurveData = prepareGalaxyForMOND(galaxyNames{ii},0.5,0.7);
end

% Create a vector with all values of a0 we want to try out:
a0Values = transpose(a0Min:a0Step:a0Max);

% Create a matrix with one cell for every combination of a value of a0 and
% a galaxy:
chiSquared = zeros(length(a0Values), numOfGalaxies);

% Iterate over all galaxies:
for jj = 1:numOfGalaxies

    % Iterate over all values of a0:
    for ii = 1:length(a0Values)

        % Calculating chi squared for the current combination of a0 and a
        % galaxy:
        chiSquared(ii,jj) = getChiSquaredForGalaxy(galaxies{jj}.rotationCurveData, a0Values(ii), interpolationFunction);
    end

    % Each galaxy gets a vector with chi squared for every value of a0:
    galaxies{jj}.chiSquared = chiSquared(:,jj);

    % Get the minimum chi squared for this galaxy:
    galaxies{jj}.chiSquaredMin = min(chiSquared(:,jj));

    % Find the index of the minimum chi squared in the vector
    % "chiSquared":
    minPos = find(galaxies{jj}.chiSquared == galaxies{jj}.chiSquaredMin);

    % Find the value of a0, that corresponds to the minimum value of chi
    % squared; this value of a0 is the best a0 value:
    galaxies{jj}.bestA0 = a0Values(minPos(1));

    % Print galaxy data to console:
    %fprintf('%s: a_0 = %d km/s^2; chi^2 = %d km^2/s^2\n', galaxies{jj}.name, galaxies{jj}.bestA0, galaxies{jj}.chiSquaredMin);
end

end


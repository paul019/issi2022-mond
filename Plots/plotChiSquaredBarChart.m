function plotChiSquaredBarChart(interpolationFunctionIds)

% Define constants:
a0Min = 0.5e-13;    % in km/s^2
a0Step = 0.1e-13;  % in km/s^2
a0Max = 2.5e-13;    % in km/s^2

a0 = 1.2e-13;       % in km/s^2

% Get the number of interpolation functions:
numOfIntFcts = length(interpolationFunctionIds);

% Initialize variables:
interpolationFunctions = cell(numOfIntFcts,1);
legendArray = cell(numOfIntFcts,1);

% Interate through all interpolation functions:
for ii = 1:numOfIntFcts
    id = interpolationFunctionIds{ii};
    interpolationFunctions{ii}.id = id;
    interpolationFunctions{ii}.name = getInterpolationFunctionName(id);
    legendArray{ii} = getInterpolationFunctionName(id);
    [galaxyNames,interpolationFunctions{ii}.galaxyFittingData] = evaluateGalaxies(a0Min,a0Step,a0Max,id);
end

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

% Initialize variables:
chiSquared = zeros(numOfGalaxies, numOfIntFcts);

% Interate through all interpolation functions:
for ii = 1:numOfIntFcts
    a0Pos = find(interpolationFunctions{ii}.galaxyFittingData{1}.a0Values == a0);

    % Iterate through all galaxies:
    for jj = 1:numOfGalaxies
        chiSquared(jj,ii) = interpolationFunctions{ii}.galaxyFittingData{jj}.chiSquared(a0Pos);
    end
end

figure
bar(chiSquared);

galaxyNames = categorical(galaxyNames);
legend(legendArray);
set(gca,'YScale','log');

end


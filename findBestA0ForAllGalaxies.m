function [a0Values, chiSquared, galaxies] = findBestA0ForAllGalaxies(a0Min,a0Step,a0Max,cleanFlag,interpolationFunction)

% Default value for interpolationFunction:
if nargin < 5
    interpolationFunction = 'linear';
end

% Import the names and metadata of all galaxies from the SPARC database:
[galaxyNames,galaxyData] = ReadLelliC;

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

% Calculate the mean chi squared for every value of a0 by averaging over
% all galaxies:
chiSquaredMean = mean(chiSquared,2);

% Get the minimum mean chi squared:
chiSquaredMeanMin = min(chiSquaredMean);

% Find the index of the minimum mean chi squared in the vector
% "chiSquaredMean":
minPos = find(chiSquaredMean == chiSquaredMeanMin);

% Find the best value of a0 for all galaxies:
bestA0Mean = a0Values(minPos(1));

% Print general findings to the console:
%fprintf('a_0 = %d; chi^2 = %d km^2/s^2\n', bestA0Mean, chiSquaredMeanMin);

%--------------------------------------------------------------------------
% Prepare figure 1:

figure('NumberTitle', 'off', 'Name', 'MOND fit');

% Create a vector with the hubble type of all galaxies:
hubbleType = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    hubbleType(jj) = galaxies{jj}.hubbleType;
end

% Create a vector with the best a0 value of all galaxies:
bestA0 = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    bestA0(jj) = galaxies{jj}.bestA0;
end

% Create a vector with the minimum chi squared for every galaxy:
chiSquaredMin = transpose(min(chiSquared));

%--------------------------------------------------------------------------
% Subplot 1 (MSWD (\chi^2) vs a_0)

if ~cleanFlag
    s = subplot(1,3,1);
else
    s = subplot(1,1,1);
end
plot(a0Values, chiSquaredMean);

set(gca,'FontSize',15);
xlabel 'a_0 [km/s^2]';
ylabel '\chi^2';
title('MSWD (\chi^2) vs a_0');
grid on;

txt = strcat('a_0 = ', sprintf('%d', bestA0Mean), ' km/s^2; \chi^2 = ', sprintf('%d', chiSquaredMeanMin));
annotation('textbox','String',txt,'Position',s.Position,'Vert','top','FitBoxToText','on','BackgroundColor','w');

% Draw one curve for each inidividual galaxy:
if ~cleanFlag
    hold on;

    for jj = 1:numOfGalaxies
        plot(a0Values, chiSquared(:,jj));
    end
end

%--------------------------------------------------------------------------
% Subplot 2 (Best value of a_0 vs hubble type)

if ~cleanFlag
    subplot(1,3,2);
    scatter(hubbleType,bestA0);
    %semilogx(chiSquaredMin,bestA0,'o');
    
    set(gca,'FontSize',15);
    xlabel 'Hubble type';
    ylabel 'a_0 [km/s^2]';
    title('Best value of a_0 vs hubble type');
    grid on;
end

%--------------------------------------------------------------------------
% Subplot 3 (MSWD (\chi^2) vs hubble type)

if ~cleanFlag
    subplot(1,3,3);
    semilogy(hubbleType,chiSquaredMin,'o');
    
    hold on;
    semilogy(hubbleType,transpose(chiSquared(minPos,:)),'o');
    
    set(gca,'FontSize',15);
    xlabel 'Hubble type';
    ylabel '\chi^2';
    title('MSWD (\chi^2) vs hubble type');
    legend('Min. \chi^2', strcat('\chi^2 for a_0 = ', sprintf('%d km/s^2', bestA0Mean)), 'Location','SouthEast')
    grid on;
end

%--------------------------------------------------------------------------

end


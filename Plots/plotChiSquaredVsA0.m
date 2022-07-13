function plotChiSquaredVsA0(galaxyNames, galaxyFittingData, cleanFlag)

[~,galaxyMetadata] = getGalaxyMetadata(galaxyNames);

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

%--------------------------------------------------------------------------
% Prepare figure 1:

figure('NumberTitle', 'off', 'Name', 'MOND fit');

% Create a vector with the hubble type of all galaxies:
hubbleType = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    hubbleType(jj) = galaxyMetadata{jj}.hubbleType;
end

% Create a vector with the best a0 value of all galaxies:
bestA0 = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    bestA0(jj) = galaxyFittingData{jj}.bestA0;
end

% Create a vector with the minimum chi squared for every galaxy:
chiSquaredMin = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    chiSquaredMin(jj) = galaxyFittingData{jj}.chiSquaredMin;
end

% Create a vector with the chi squared of every galaxy for the best a0
% overall:
chiSquaredForBestOverallA0 = zeros(numOfGalaxies,1);
bestOverallA0_index = galaxyFittingData{end}.bestA0_index;
for jj = 1:numOfGalaxies
    chiSquaredForBestOverallA0(jj) = galaxyFittingData{jj}.chiSquared(bestOverallA0_index);
end

%--------------------------------------------------------------------------
% Subplot 1 (MSWD (\chi^2) vs a_0)

if ~cleanFlag
    s = subplot(1,3,1);
else
    s = subplot(1,1,1);
end
plot(galaxyFittingData{end}.a0Values, galaxyFittingData{end}.chiSquared);

set(gca,'FontSize',15);
xlabel 'a_0 [km/s^2]';
ylabel '\chi^2';
title('MSWD (\chi^2) vs a_0');
grid on;

txt = strcat('a_0 = ', sprintf('%d', galaxyFittingData{end}.bestA0), ' km/s^2; \chi^2 = ', sprintf('%d', galaxyFittingData{end}.chiSquaredMin));
annotation('textbox','String',txt,'Position',s.Position,'Vert','top','FitBoxToText','on','BackgroundColor','w');

% Draw one curve for each inidividual galaxy:
if ~cleanFlag
    hold on;

    for jj = 1:numOfGalaxies
        plot(a0Values, galaxyFittingData{jj}.chiSquared);
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
    semilogy(hubbleType,chiSquaredForBestOverallA0,'o');
    
    set(gca,'FontSize',15);
    xlabel 'Hubble type';
    ylabel '\chi^2';
    title('MSWD (\chi^2) vs hubble type');
    legend('Min. \chi^2', strcat('\chi^2 for a_0 = ', sprintf('%d km/s^2', bestA0Mean)), 'Location','SouthEast')
    grid on;
end

%--------------------------------------------------------------------------

end


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
chiSquaredReducedMin = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    chiSquaredReducedMin(jj) = galaxyFittingData{jj}.chiSquaredReducedMin;
end

% Create a vector with the chi squared of every galaxy for the best a0
% overall:
chiSquaredReducedForBestOverallA0 = zeros(numOfGalaxies,1);
bestOverallA0_index = galaxyFittingData{end}.bestA0_index;
for jj = 1:numOfGalaxies
    chiSquaredReducedForBestOverallA0(jj) = galaxyFittingData{jj}.chiSquaredReduced(bestOverallA0_index);
end

%--------------------------------------------------------------------------
% Subplot 1 (MSWD (\chi_\nu^2) vs a_0)

if ~cleanFlag
    s = subplot(1,3,1);
else
    s = subplot(1,1,1);
end
plot(galaxyFittingData{end}.a0Values, galaxyFittingData{end}.chiSquaredReduced);

set(gca,'FontSize',15);
xlabel 'a_0 [km/s^2]';
ylabel '\chi_\nu^2';
title('MSWD (\chi_\nu^2) vs a_0');
grid on;

txt = strcat('a_0 = ', sprintf('%d', galaxyFittingData{end}.bestA0), ' km/s^2; \chi_\nu^2 = ', sprintf('%d', galaxyFittingData{end}.chiSquaredReducedMin));
annotation('textbox','String',txt,'Position',s.Position,'Vert','top','FitBoxToText','on','BackgroundColor','w');

% Draw one curve for each inidividual galaxy:
if ~cleanFlag
    hold on;

    for jj = 1:numOfGalaxies
        plot(galaxyFittingData{jj}.a0Values, galaxyFittingData{jj}.chiSquaredReduced);
    end
end

%--------------------------------------------------------------------------
% Subplot 2 (Best value of a_0 vs hubble type)

if ~cleanFlag
    subplot(2,3,2);
    scatter(hubbleType,bestA0);
    
    set(gca,'FontSize',15);
    xlabel 'Hubble type';
    ylabel 'Best value of a_0 [km/s^2]';
    title('Best value of a_0 vs hubble type');
    grid on;
end

%--------------------------------------------------------------------------
% Subplot 3 (MSWD (\chi_\nu^2) vs hubble type)

if ~cleanFlag
    subplot(2,3,3);
    semilogy(hubbleType,chiSquaredReducedMin,'o');
    
    hold on;
    semilogy(hubbleType,chiSquaredReducedForBestOverallA0,'o');
    
    set(gca,'FontSize',15);
    xlabel 'Hubble type';
    ylabel '\chi_\nu^2';
    title('MSWD (\chi_\nu^2) vs hubble type');
    legend('Min. \chi_\nu^2', strcat('\chi_\nu^2 for a_0 = ', sprintf('%d km/s^2', galaxyFittingData{end}.bestA0)), 'Location','SouthEast')
    grid on;
end

%--------------------------------------------------------------------------
% Subplot 4 (MSWD (\chi_\nu^2) vs best value of a0)

if ~cleanFlag
    subplot(2,3,[5,6]);
    semilogy(bestA0,chiSquaredReducedMin,'o');
    
    hold on;
    semilogy(bestA0,chiSquaredReducedForBestOverallA0,'o');
    
    set(gca,'FontSize',15);
    xlabel 'Best value of a_0 [km/s^2]';
    ylabel '\chi_\nu^2';
    title('MSWD (\chi_\nu^2) vs best value of a_0');
    legend('\chi_\nu^2 for respective a_0', strcat('\chi_\nu^2 for a_0 = ', sprintf('%d km/s^2', galaxyFittingData{end}.bestA0)), 'Location','SouthEast')
    grid on;
end

%--------------------------------------------------------------------------

end


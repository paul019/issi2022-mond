function plotChiSquaredVsA0(a0Values, chiSquared, galaxies, cleanFlag)

% Get the number of galaxies:
numOfGalaxies = length(galaxies);

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


function [a0Values, chiSquaredVector] = findBestA0ForAllGalaxies(a0Min,a0Step,a0Max,myGalaxyNames)

[galaxyNames,galaxyData_] = ReadLelliC;

if exist('galaxyNames','var')
    myGalaxyNames = galaxyNames;
    galaxyData = galaxyData_;
else
    for ii=1:length(galaxyNames)
        if galaxyNames{ii} 

        end
    end
end

numOfGalaxies = length(galaxyNames);
galaxies = cell(numOfGalaxies,1);

for ii = 1:numOfGalaxies
    galaxies{ii}.name = galaxyNames{ii};
    galaxies{ii}.hubbleType = galaxyData{ii}(1);
    galaxies{ii}.rotationCurveData = prepareGalaxyForMOND(galaxyNames{ii},0.5,0.5);
end

a0Values = transpose(a0Min:a0Step:a0Max);
chiSquared = zeros(length(a0Values), numOfGalaxies);

iterationCounter = 0;
numOfIterations = numOfGalaxies * length(a0Values);

for jj = 1:numOfGalaxies
    rotationCurveData = galaxies{jj}.rotationCurveData;

    for ii = 1:length(a0Values)
        chiSquared(ii,jj) = getChiSquaredForGalaxy(rotationCurveData, a0Values(ii));

        iterationCounter = iterationCounter + 1;
        %progmeter(iterationCounter/numOfIterations, 'Finding the best value for a0');
    end

    galaxies{jj}.chiSquared = chiSquared(:,jj);
    galaxies{jj}.chiSquaredMin = min(chiSquared(:,jj));
    minPos = find(chiSquared(:,jj) == galaxies{jj}.chiSquaredMin);
    galaxies{jj}.bestA0 = a0Values(minPos(1));

    %fprintf('%s: a_0 = %d km/s^2; chi^2 = %d km^2/s^2\n', galaxies{jj}.name, galaxies{jj}.bestA0, galaxies{jj}.chiSquaredMin);
end

chiSquaredMean = transpose(mean(transpose(chiSquared)));
chiSquaredMeanMin = min(chiSquaredMean);
minPos = find(chiSquaredMean == chiSquaredMeanMin);
bestA0Mean = a0Values(minPos(1));

%fprintf('a_0 = %d; chi^2 = %d km^2/s^2\n', bestA0Mean, chiSquaredMeanMin);

%--------------------------------------------------------------------------

figure('NumberTitle', 'off', 'Name', 'MOND fit');

hubbleType = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    hubbleType(jj) = galaxies{jj}.hubbleType;
end
bestA0 = zeros(numOfGalaxies,1);
for jj = 1:numOfGalaxies
    bestA0(jj) = galaxies{jj}.bestA0;
end
chiSquaredMin = transpose(min(chiSquared));

%--------------------------------------------------------------------------

s = subplot(1,3,1);

plot(a0Values, chiSquaredMean);

set(gca,'FontSize',15);
xlabel 'a_0 [km/s^2]';
ylabel '\chi^2 [km^2/s^2]';
title('MSWD (\chi^2) vs a_0');
grid on;

txt = strcat('a_0 = ', sprintf('%d', bestA0Mean), ' km/s^2; \chi^2 = ', sprintf('%d', chiSquaredMeanMin), ' km^2/s^2');
annotation('textbox','String',txt,'Position',s.Position,'Vert','top','FitBoxToText','on','BackgroundColor','w')

%--------------------------------------------------------------------------

subplot(1,3,2);
scatter(hubbleType,bestA0);

set(gca,'FontSize',15);
xlabel 'Hubble type';
ylabel 'a_0 [km/s^2]';
title('Best value of a_0 vs hubble type');
grid on;

%--------------------------------------------------------------------------

subplot(1,3,3);
semilogy(hubbleType,chiSquaredMin,'o');

hold on;
semilogy(hubbleType,transpose(chiSquared(minPos,:)),'o');

set(gca,'FontSize',15);
xlabel 'Hubble type';
ylabel '\chi^2 [km^2/s^2]';
title('MSWD (\chi^2) vs hubble type');
legend('Min. \chi^2', strcat('\chi^2 for a_0 = ', sprintf('%d km/s^2', bestA0Mean)), 'Location','SouthEast')
grid on;

%--------------------------------------------------------------------------

%{
hold on;

for jj = 1:numOfGalaxies
    plot(a0Values, chiSquared(:,jj));
end

set(gca,'FontSize',15);
xlabel 'a_0 [km/s^2]';
ylabel '\chi^2 [km^2/s^2]';

%}

end


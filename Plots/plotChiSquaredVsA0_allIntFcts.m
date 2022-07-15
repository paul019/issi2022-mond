function plotChiSquaredVsA0_allIntFcts(a0Min,a0Step,a0Max)

interpolationFunctions = {
    struct('id','linear'),...
    struct('id','rar'),...
    struct('id','simple'),...
    struct('id','standard'),...
    struct('id','toy'),...
    struct('id','exp')
};

for ii = 1:length(interpolationFunctions)
    id = interpolationFunctions{ii}.id;
    interpolationFunctions{ii}.name = getInterpolationFunctionName(id);
    [~,interpolationFunctions{ii}.galaxyFittingData] = evaluateGalaxies(a0Min,a0Step,a0Max,id);
end

figure('NumberTitle', 'off', 'Name', 'MOND fit');

legendText = cell(length(interpolationFunctions),1);

for ii = 1:length(interpolationFunctions)
    galaxyFittingData = interpolationFunctions{ii}.galaxyFittingData;

    plot(galaxyFittingData{end}.a0Values, galaxyFittingData{end}.chiSquaredReduced);
    hold on;

    bestA0 = galaxyFittingData{end}.bestA0;
    chiSquaredMin = galaxyFittingData{end}.chiSquaredMin;

    legendText{ii} = strcat(interpolationFunctions{ii}.name, ' (a_0 = ', num2str(bestA0), ' ms^{-2})');
end

% Get the layout right:
set(gca,'FontSize',15);
xlabel 'a_0 [km/s^2]';
ylabel '\chi_\nu^2';
title('MSWD (\chi_\nu^2) vs a_0');
grid on;
legend(legendText, 'Location', 'NorthEast');

end


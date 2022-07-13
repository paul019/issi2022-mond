function plotChiSquaredVsA0_allIntFcts(a0Min,a0Step,a0Max)

interpolationFunctions = {
    struct('id','linear','name','Linear approach'),...
    struct('id','rar','name','Radial acceleration relation (RAR)'),...
    struct('id','simple','name','Simple interpolation'),...
    struct('id','standard','name','Standard interpolation'),...
    struct('id','toy','name','Toy'),...
    struct('id','exp','name','''Exponential'' interpolation')
};

for ii = 1:length(interpolationFunctions)
    id = interpolationFunctions{ii}.id;
    [~,interpolationFunctions{ii}.galaxyFittingData] = evaluateGalaxies(a0Min,a0Step,a0Max,id);
end

figure('NumberTitle', 'off', 'Name', 'MOND fit');

legendText = cell(length(interpolationFunctions),1);

for ii = 1:length(interpolationFunctions)
    galaxyFittingData = interpolationFunctions{ii}.galaxyFittingData;

    plot(galaxyFittingData{end}.a0Values, galaxyFittingData{end}.chiSquared);
    hold on;

    legendText{ii} = interpolationFunctions{ii}.name;
end

% Get the layout right:
set(gca,'FontSize',15);
xlabel 'a_0 [km/s^2]';
ylabel '\chi^2';
title('MSWD (\chi^2) vs a_0');
grid on;
legend(legendText, 'Location', 'NorthEast');

end


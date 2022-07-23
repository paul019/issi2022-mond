function plotChiSquaredVsA0_multipleIntFcts(intFctIds,galaxyFittingDataArray)

figure('NumberTitle', 'off', 'Name', 'MOND fit');

legendText = cell(length(intFctIds),1);

for ii = 1:length(intFctIds)
    galaxyFittingData = galaxyFittingDataArray{ii};

    plot(galaxyFittingData{end}.a0Values * 10^3, galaxyFittingData{end}.chiSquaredReduced);
    hold on;

    bestA0 = galaxyFittingData{end}.bestA0;
    chiSquaredReducedMin = galaxyFittingData{end}.chiSquaredReducedMin;

    legendText{ii} = [getInterpolationFunctionName(intFctIds{ii}), ' (a_0 = ', num2str(bestA0 * 1e3), ' m/s^2)'];
end

% Get the layout right:
set(gca,'FontSize',15);
xlabel 'a_0 [m/s^2]';
ylabel '\chi_\nu^2';
ylim([0 200])
title('MSWD per degree of freedom (\chi_\nu^2) vs a_0');
grid on;
legend(legendText, 'Location', 'NorthEast');

end


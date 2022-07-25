function plotChiSquaredVsA0_multipleIntFcts(intFctIds,galaxyFittingDataArray,easyUnitsFlag,yRange)

if nargin < 3
    easyUnitsFlag = false;
end

figure;

legendText = cell(length(intFctIds),1);

for ii = 1:length(intFctIds)
    galaxyFittingData = galaxyFittingDataArray{ii};

    if easyUnitsFlag
        galaxyFittingData{end}.a0Values = galaxyFittingData{end}.a0Values * 10^12;  % convert to nm/s^2
        galaxyFittingData{end}.bestA0   = galaxyFittingData{end}.bestA0 * 10^12;    % convert to nm/s^2
    else
        galaxyFittingData{end}.a0Values = galaxyFittingData{end}.a0Values * 10^3;   % convert to m/s^2
        galaxyFittingData{end}.bestA0   = galaxyFittingData{end}.bestA0 * 10^3;     % convert to m/s^2
    end

    p = plot(galaxyFittingData{end}.a0Values, galaxyFittingData{end}.chiSquaredReduced);
    p.Color = getInterpolationFunctionColor(intFctIds{ii});
    p.LineWidth = 2;
    hold on;

    bestA0 = galaxyFittingData{end}.bestA0;
    %chiSquaredReducedMin = galaxyFittingData{end}.chiSquaredReducedMin;

    if easyUnitsFlag
        legendText{ii} = [getInterpolationFunctionName(intFctIds{ii}), ' (a_0 = ', num2str(bestA0), ' nm/s^2)'];
    else
        legendText{ii} = [getInterpolationFunctionName(intFctIds{ii}), ' (a_0 = ', num2str(bestA0), ' m/s^2)'];
    end
end

% Get the layout right:
set(gca,'FontSize',15);
if easyUnitsFlag
    xlabel 'a_0 [nm/s^2]';
else
    xlabel 'a_0 [m/s^2]';
end
ylabel '\chi_\nu^2';
if nargin >= 4
    ylim(yRange)
end
title('MSWD per degree of freedom (\chi_\nu^2) vs a_0');
grid on;
legend(legendText, 'Location', 'NorthWest');

end


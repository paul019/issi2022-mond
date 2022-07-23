function [fitIds,galaxyFittingDataArray] = exportFits

%--------------------------------------------------------------------------
% Preparations:
%--------------------------------------------------------------------------

fprintf('\n');

% Define constants:
a0Step = 0.1e-13;       % in km/s^2
a0Min  = a0Step;         % in km/s^2
a0Max  = 2.5e-13;        % in km/s^2
minQuality = 3;

% Array of galaxies to be used:
galaxyNames = getGalaxiesByQuality(minQuality);
numOfGalaxies = length(galaxyNames);

% Array of galaxies to be plotted:
galaxyNamesForPlot = {'UGC01281','NGC0024','NGC5055'};

% Array of interpolation functions and fits to be used:
%intFctIds = { 'linear'; 'rar'; 'simple'; 'standard'; 'toy'; 'exp' };
intFctIds = { 'linear'; 'rar'; 'simple'; 'standard' };
fitIds = [intFctIds; {'nfw'}];

numOfIntFcts = length(intFctIds);
numOfFits = numOfIntFcts + 1;

fitNames = cell(numOfFits, 1);
for ii = 1:numOfIntFcts
    fitNames{ii} = getInterpolationFunctionName(intFctIds{ii});
end
fitNames{end} = 'NFW fit';

% Initialize array of galaxy fitting data for different interpolation
% functions and NFW fit:
galaxyFittingDataArray = cell(numOfFits, 1);

%--------------------------------------------------------------------------
% MOND fit:
%--------------------------------------------------------------------------

% Iterate through all interpolation functions:
for ii = 1:numOfIntFcts
    fprintf('Evaluating interpolation function %d of %d\n', ii, length(intFctIds))

    [~, galaxyFittingDataArray{ii}] = evaluateGalaxiesForMOND(a0Min, a0Step, a0Max, intFctIds{ii}, false, galaxyNames);
end

%--------------------------------------------------------------------------
% NFW fit:
%--------------------------------------------------------------------------

fprintf('Evaluating NFW fit\n');

[~, galaxyFittingDataArray{end}] = evaluateGalaxiesForNFW(false, galaxyNames);

%--------------------------------------------------------------------------
% Status update:
%--------------------------------------------------------------------------

fprintf('\nEvaluated %d datapoints from %d galaxies.\n', galaxyFittingDataArray{ii}{end}.numberOfDatapoints, numOfGalaxies);

fprintf('Exporting results\n')

%--------------------------------------------------------------------------
% Export first table (MONDFits_overview.csv):
%--------------------------------------------------------------------------

% Initialize more arrays to save data:
bestA0 = zeros(numOfFits, 1);
chiSquaredReduced = zeros(numOfFits, numOfGalaxies + 1);
numOfGalaxiesWithGoodFit = zeros(numOfFits, 1);

% Iterate through all interpolation functions:
for ii = 1:numOfIntFcts
    bestA0(ii) = galaxyFittingDataArray{ii}{end}.bestA0 * 1e3;
end

% Iterate through all fits:
for ii = 1:numOfFits
    for jj = 1:numOfGalaxies+1
        chiSquaredReduced(ii,jj) = galaxyFittingDataArray{ii}{jj}.chiSquaredReduced_general;
    end

    % Iterate thorugh all galaxies:
    for jj = 1:numOfGalaxies
        if chiSquaredReduced(ii,jj) < 1
            numOfGalaxiesWithGoodFit(ii) = numOfGalaxiesWithGoodFit(ii) + 1;
        end
    end
end

column1 = fitNames;
column2 = bestA0;
column3 = chiSquaredReduced(:,numOfGalaxies + 1);
column4 = numOfGalaxiesWithGoodFit;

% Create and save table:
table1 = table(column1, column2, column3, column4, 'VariableNames', {'Fit', 'Best a_0 (in m/s^2)', 'chi_v^2', 'Number of galaxies with chi_v^2 < 1'});

writetable(table1,'output/MONDFits_overview.csv');

%--------------------------------------------------------------------------
% Export second table (MONDFits_galaxies.csv):
%--------------------------------------------------------------------------

% Initialize row array with interpolation function names:
fitNames_ = cell(1, numOfFits);
for ii = 1:numOfFits
    fitNames_{ii} = ['chi_v^2 for ', fitNames{ii}];
end

% Initialize column arrays with data:
galaxyNames_ = [galaxyNames; {'All galaxies'}];
numberOfDatapoints = zeros(numOfGalaxies + 1, 1);
chiSquaredReduced_ = num2cell(chiSquaredReduced',1);
p0 = zeros(numOfGalaxies + 1, 1);
R_S = zeros(numOfGalaxies + 1, 1);

% Iterate through all galaxies:
for ii = 1:numOfGalaxies
    numberOfDatapoints(ii) = galaxyFittingDataArray{1}{ii}.numberOfDatapoints;

    p0(ii) = galaxyFittingDataArray{end}{ii}.bestP0;
    R_S(ii) = galaxyFittingDataArray{end}{ii}.bestR_S;
end
numberOfDatapoints(end) = galaxyFittingDataArray{1}{end}.numberOfDatapoints;

% Sort the galaxies by chi squared reduced for the simple fit:
[~, order] = sort(chiSquaredReduced(3,1:numOfGalaxies));
order(numOfGalaxies + 1) = numOfGalaxies + 1;
column1 = galaxyNames_(order);
column2 = numberOfDatapoints(order);
columns3 = cell(numOfFits,1);
for ii = 1:numOfFits
    columns3{ii} = chiSquaredReduced_{ii}(order);
end
column4 = p0(order);
column5 = R_S(order);

% Create and save table:

table2 = table(column1, column2, columns3{:}, column4, column5, 'VariableNames', [{'Galaxy Name', 'Number of datapoints'}, fitNames_, {'\rho_0', 'R_S'}]);

writetable(table2,'output/MONDFits_galaxies.csv');

%--------------------------------------------------------------------------
% Export first plot (MSWD_vs_a0.png)
%--------------------------------------------------------------------------

plotChiSquaredVsA0_multipleIntFcts(intFctIds,galaxyFittingDataArray);

set(gcf, 'Position',  [100, 100, 1000, 750])
saveas(gcf,'output/MSWD_vs_a0.png');

%--------------------------------------------------------------------------
% Export second plot (MSWD_histogram.png)
%--------------------------------------------------------------------------

figure

for ii = 1:numOfFits
    subplot(2,4,ii);

    title = fitNames{ii};
    if ii ~= numOfFits
        subtitle = sprintf('a_0 = %s m/s^2', num2str(galaxyFittingDataArray{ii}{end}.bestA0 * 1e3));
    else
        subtitle = '';
    end

    plotChiSquaredHistogram(galaxyNames,galaxyFittingDataArray{ii},50,false,title,subtitle);
end

set(gcf, 'Position',  [100, 100, 1500, 750])
saveas(gcf,'output/MONDFits_MSWD_histogram.png');

%--------------------------------------------------------------------------
% Export plots for single galaxies (GALAXYNAME_v_vs_r.png)
%--------------------------------------------------------------------------

for galaxyName_ = galaxyNamesForPlot
    galaxyName = galaxyName_{1};
    galaxyIndex = find(strcmp(galaxyNames, galaxyName));

    % Prepare MOND fits:
    mondFits = cell(numOfIntFcts,1);
    for ii = 1:numOfIntFcts
        mondFits{ii} = struct('intFctId',intFctIds{ii},'a0',galaxyFittingDataArray{ii}{end}.bestA0,'chiSquaredReduced',galaxyFittingDataArray{ii}{galaxyIndex}.chiSquaredReduced_general);
    end

    % Prepare NFW fit:
    nfwFits = {struct('p0',p0(galaxyIndex),'R_S',R_S(galaxyIndex),'chiSquaredReduced',galaxyFittingDataArray{end}{galaxyIndex}.chiSquaredReduced_general)};

    plotGalaxyVelocityWithFits(galaxyName, 0.5, 0.7, mondFits, nfwFits, false, false);

    set(gcf, 'Position',  [100, 100, 1000, 750])
    saveas(gcf,sprintf('output/%s_v_vs_r.png', galaxyName));
end

fprintf('\n');

end


function [fitIds,galaxyFittingDataArray] = exportFits(galaxyNamesForPlot, intFctIds, minQuality)

%--------------------------------------------------------------------------
% Preparations:
%--------------------------------------------------------------------------

fprintf('\n');

% Create output folder:
if not(isfolder('output'))
   mkdir('output')
end

% Define constants:
a0Step = 0.01e-13;       % in km/s^2
a0Min  = a0Step;         % in km/s^2
a0Max  = 2.5e-13;        % in km/s^2
if nargin < 3
    minQuality = 3;
end

% Array of galaxies to be used:
galaxyNames = getGalaxiesByQuality(minQuality);
[~, galaxyMetadata] = getGalaxyMetadata(galaxyNames);
numOfGalaxies = length(galaxyNames);

% Array of galaxies to be plotted:
if nargin < 1
    galaxyNamesForPlot = {'UGC01281','NGC6503','NGC5055'};
end

% Array of interpolation functions and fits to be used:
if nargin < 2
    intFctIds = getAllInterpolationFunctionIds;
end
%intFctIds = { 'linear'; 'rar'; 'simple'; 'standard' };
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
p0_eV = zeros(numOfGalaxies + 1, 1);
R_S = zeros(numOfGalaxies + 1, 1);
R_S_kpc = zeros(numOfGalaxies + 1, 1);

% Iterate through all galaxies:
for ii = 1:numOfGalaxies
    numberOfDatapoints(ii) = galaxyFittingDataArray{1}{ii}.numberOfDatapoints;

    p0(ii) = galaxyFittingDataArray{end}{ii}.bestP0;
    p0_eV(ii) = galaxyFittingDataArray{end}{ii}.bestP0_eV;
    R_S(ii) = galaxyFittingDataArray{end}{ii}.bestR_S;
    R_S_kpc(ii) = galaxyFittingDataArray{end}{ii}.bestR_S_kpc;
end
numberOfDatapoints(end) = galaxyFittingDataArray{1}{end}.numberOfDatapoints;

% Sort the galaxies by chi squared reduced for the simple fit if possible:
[~, order] = sort(chiSquaredReduced(min(3,numOfFits),1:numOfGalaxies));
order(numOfGalaxies + 1) = numOfGalaxies + 1;
column1 = galaxyNames_(order);
column2 = numberOfDatapoints(order);
column3 = [cellfun(@(x) num2str(x.qualityFlag), galaxyMetadata); '-'];
column3 = column3(order);
columns4 = cell(numOfFits,1);
for ii = 1:numOfFits
    columns4{ii} = chiSquaredReduced_{ii}(order);
end
column5 = p0_eV(order);
column5 = arrayfun(@(x) num2str(x), column5,'UniformOutput',false);
column5{end} = '-';
column6 = R_S_kpc(order);
column6 = arrayfun(@(x) num2str(x), column6,'UniformOutput',false);
column6{end} = '-';

disp(size(column1));
disp(size(column2));
disp(size(column3));
disp(size(columns4{2}));
disp(size(column5));
disp(size(column6));

% Create and save table:

table2 = table(column1, column2, column3, columns4{:}, column5, column6, 'VariableNames', [{'Galaxy Name', 'Number of datapoints', 'Quality Flag'}, fitNames_, {'\rho_0 [eV/cm^3]', 'R_S [kpc]'}]);

writetable(table2,'output/MONDFits_galaxies.csv');

%--------------------------------------------------------------------------
% Export first plot (MSWD_vs_a0.png)
%--------------------------------------------------------------------------

plotChiSquaredVsA0_multipleIntFcts(intFctIds,galaxyFittingDataArray,false,[40,120]);

set(gcf, 'Position',  [100, 100, 1000, 500])
saveas(gcf,'output/MSWD_vs_a0.png');

%--------------------------------------------------------------------------
% Export second plot (MSWD_histogram.png)
%--------------------------------------------------------------------------

figure

for ii = 1:numOfFits
    subplot(2,4,ii);

    title_ = fitNames{ii};
    if ii ~= numOfFits
        subtitle = sprintf('a_0 = %s m/s^2', num2str(galaxyFittingDataArray{ii}{end}.bestA0 * 1e3));
    else
        subtitle = '';
    end

    plotChiSquaredHistogram(galaxyNames,galaxyFittingDataArray{ii},50,false,title_,subtitle);
end

set(gcf, 'Position',  [100, 100, 1500, 750])
saveas(gcf,'output/MONDFits_MSWD_histogram.png');

%--------------------------------------------------------------------------
% Export third plot (intFcts.png)
%--------------------------------------------------------------------------

plotInterpolationFunctions;

set(gcf, 'Position',  [100, 100, 500, 400])
saveas(gcf,'output/intFcts.png');

%--------------------------------------------------------------------------
% Export plots for single galaxies (GALAXYNAME_v_vs_r.png)
%--------------------------------------------------------------------------

for galaxyName_ = galaxyNamesForPlot
    galaxyName = galaxyName_{1};
    galaxyIndex = find(strcmp(galaxyNames, galaxyName));

    % Prepare MOND fits:
    mondFits = cell(numOfIntFcts,1);
    for ii = 1:numOfIntFcts
        mondFits{ii} = galaxyFittingDataArray{ii}{galaxyIndex}.mondFit;
    end

    % Prepare NFW fit:
    nfwFits = {galaxyFittingDataArray{end}{galaxyIndex}.nfwFit};

    plotGalaxyRotationCurveWithFits(galaxyName, 0.5, 0.7, mondFits, nfwFits, false, false);

    set(gcf, 'Position',  [100, 100, 1000, 750])
    saveas(gcf,sprintf('output/%s_v_vs_r.png', galaxyName));
end

fprintf('\n');

end


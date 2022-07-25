function [galaxyNames, galaxyFittingData] = evaluateGalaxiesForNFW(printflag,galaxyNames)

% Define constants:
MtoLdisk = 0.5;
MtoLbulge = 0.7;
G = 6.6743e-2;                      % in km^3 * kg^(-1) * s^(-2)
c = 299792458;                      % in m/s
jouleInEV = 6.242e18;
kpcInKm = 3.0857e16;

% Define starting values for the fit with the right order of magnitude:
params_start = [1e-23, 1e17];
    % params_start(1) = p0 in kg/km^3
    % params_start(2) = R_S in km

% Import the names and metadata of all galaxies from the SPARC database:
[allGalaxyNames,~] = ReadLelliC;

% Deal with the input argument "galaxyNames":
if nargin < 2
    galaxyNames = allGalaxyNames;
else
    myGalaxyNames = galaxyNames;

    galaxyNames = {};
    index = 1;

    for ii = 1:length(allGalaxyNames)
        for jj = 1:length(myGalaxyNames)
            if strcmp(allGalaxyNames{ii}, myGalaxyNames{jj})
                galaxyNames{index} = allGalaxyNames{ii};
                index = index + 1;
                break
            end
        end
    end
end

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

% Create a cell array for all galaxies:
galaxyFittingData = cell(numOfGalaxies + 1, 1);

% Keep track of the total number of datapoints:
totalNumberOfDatapoints = 0;

% Initialize a vector for all chi squared values:
chiSquared = zeros(numOfGalaxies,1);

%--------------------------------------------------------------------------

% Iterate through all galaxies in order to get rotation curve data:
for ii = 1:numOfGalaxies
    galaxyFittingData{ii}.rotationCurveData = prepareGalaxyRotationCurveData(galaxyNames{ii},MtoLdisk,MtoLbulge,true);
    totalNumberOfDatapoints = totalNumberOfDatapoints + length(galaxyFittingData{ii}.rotationCurveData);

    galaxyFittingData{ii}.typeOfFit = 'NFW';
    galaxyFittingData{ii}.name = galaxyNames{jj};
end

% Print status update:
if printflag
    fprintf('\nEvaluating %d datapoints from %d galaxies.\n', totalNumberOfDatapoints, numOfGalaxies);
end

% Define fitting functions:
V_DM_fun       = @(r,params) sqrt( G * r.^(-1) * 4 * pi * params(1) * params(2)^3 .* ( log( (params(2) + r) ./ params(2) ) + params(2) ./ (params(2) + r) - 1 ) );
V_total_fun    = @(r,V_bar,params) sqrt( V_bar.^2 + V_DM_fun(r,params).^2 );

    function chiSquared = chiSquared_fun(params,data)
        if params(1) < 0 || params(2) < 0
            chiSquared = 1e10;
        else
            chiSquared = sum( ( ( V_total_fun(data(:,9),data(:,11),params) - data(:,2) ) ./ data(:,3) ).^2 );
        end
    end

%--------------------------------------------------------------------------

% Iterate over all galaxies:
for jj = 1:numOfGalaxies
    % Use the built-in MATLAB function fminsearch in order to find the best
    % parameters p0 and R_S:
    options = optimset('MaxFunEvals',1000000,'MaxIter',1000000);
    [params_best,chiSquared_min] = fminsearch(@chiSquared_fun, params_start, options, galaxyFittingData{jj}.rotationCurveData);

    %{
    if min(V_DM_fun(galaxyFittingData{jj}.rotationCurveData(:,9), params_best)) < 1
        disp(jj)
        disp(galaxyNames{jj})
        disp(V_DM_fun(galaxyFittingData{jj}.rotationCurveData(:,9), params_best))
        disp(params_best)
        disp(chiSquared_min)
    end
    %}

    % Calculate the degrees of freedom:
    galaxyFittingData{jj}.degreesOfFreedom = length(galaxyFittingData{jj}.rotationCurveData) - 3;

    % Save chi squared data:
    galaxyFittingData{jj}.chiSquared = chiSquared_min;
    chiSquared(jj) = chiSquared_min;
    galaxyFittingData{jj}.chiSquaredReduced = chiSquared_min / galaxyFittingData{jj}.degreesOfFreedom;

    galaxyFittingData{jj}.chiSquared_general = galaxyFittingData{jj}.chiSquared;
    galaxyFittingData{jj}.chiSquaredReduced_general = galaxyFittingData{jj}.chiSquaredReduced;

    % Save parameter data:
    galaxyFittingData{jj}.bestP0 = params_best(1);      % in kg/km^3
    galaxyFittingData{jj}.bestR_S = params_best(2);     % in km

    galaxyFittingData{jj}.bestP0_eV = params_best(1) / 1e9 * c^2 * jouleInEV / 1e6;       % in eV/cm^3
    galaxyFittingData{jj}.bestR_S_kpc = params_best(2) / kpcInKm;                         % in kpc

    % Save galaxy data as human-readable string:
    galaxyFittingData{jj}.dataString = sprintf('%s: p_0 = %d kg/km^3; R_S = %d km; chi_v^2 = %d', galaxyNames{jj}, galaxyFittingData{jj}.bestP0, galaxyFittingData{jj}.bestR_S, galaxyFittingData{jj}.chiSquaredReduced);

    % Save galaxy data as nfw fit (input to
    % "plotGalaxyRotationCurveWithFits"):
    galaxyFittingData{jj}.nfwFit = struct('p0',params_best(1),'R_S',params_best(2),'chiSquaredReduced',galaxyFittingData{jj}.chiSquaredReduced_general);
end

%--------------------------------------------------------------------------

% The last entry in the cell array "galaxies" represents the average of all
% galaxies:
jj = numOfGalaxies + 1;

galaxyFittingData{jj}.typeOfFit = 'MOND';
galaxyFittingData{jj}.name = 'All galaxies';

% Calculate the degrees of freedom:
galaxyFittingData{jj}.degreesOfFreedom = totalNumberOfDatapoints - 2 * numOfGalaxies - 1;

% Calculate the chi squared for all galaxies:
galaxyFittingData{jj}.chiSquared = sum(chiSquared);
galaxyFittingData{jj}.chiSquaredReduced = sum(chiSquared) / galaxyFittingData{jj}.degreesOfFreedom;

galaxyFittingData{jj}.chiSquared_general = galaxyFittingData{jj}.chiSquared;
galaxyFittingData{jj}.chiSquaredReduced_general = galaxyFittingData{jj}.chiSquaredReduced;

% Print general findings to the console:
if printflag
    fprintf('All galaxies: chi_v^2 = %d\n\n', galaxyFittingData{jj}.chiSquaredReduced);
end

%--------------------------------------------------------------------------
% Print galaxy data to console:

for jj = 1:numOfGalaxies
    if printflag
        fprintf('%s\n', galaxyFittingData{jj}.dataString);
    end
end

%--------------------------------------------------------------------------
% Determine the galaxy with the best fit:

bestChiSquaredReduced = 10^10;
bestGalaxyName = '';
bestGalaxyIndex = -1;

for jj = 1:numOfGalaxies
    galaxy_chiSquaredReduced = galaxyFittingData{jj}.chiSquaredReduced;

    if galaxy_chiSquaredReduced < bestChiSquaredReduced
        bestChiSquaredReduced = galaxy_chiSquaredReduced;
        bestGalaxyName = galaxyNames{jj};
        bestGalaxyIndex = jj;
    end
end

if printflag
    fprintf('Best galaxy: %s\n\n', galaxyFittingData{bestGalaxyIndex}.dataString);
end

end


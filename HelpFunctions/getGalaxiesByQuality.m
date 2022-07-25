function galaxyNames = getGalaxiesByQuality(minQuality)
% quality: 1 – High; 2 – Medium; 3 – Low

% Get all galaxies:
[allGalaxyNames,allGalaxyData] = ReadLelliC;

% Number of galaxies with sufficient quality:
numOfQualityGalaxies = 0;

% Iterate through all galaxies:
for ii = 1:length(allGalaxyNames)
    if allGalaxyData{ii}(17) <= minQuality
        numOfQualityGalaxies = numOfQualityGalaxies + 1;
    end
end

% Initialize galaxy names array:
galaxyNames = cell(numOfQualityGalaxies, 1);

% Iterate through all galaxies:
jj = 1;
for ii = 1:length(allGalaxyNames)
    % Only add galaxies with sufficient quality:
    if allGalaxyData{ii}(17) <= minQuality
        galaxyNames{jj} = allGalaxyNames{ii};
        jj = jj + 1;
    end
end

end


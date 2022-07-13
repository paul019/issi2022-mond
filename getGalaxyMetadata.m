function [galaxyNames,galaxyMetadata] = getGalaxyMetadata(galaxyNames)

% Import the names and metadata of all galaxies from the SPARC database:
[allGalaxyNames,allGalaxyData] = ReadLelliC;

% Deal with the input argument "galaxyNames":
if nargin < 1
    galaxyNames = allGalaxyNames;
    galaxyData = allGalaxyData;
else
    myGalaxyNames = galaxyNames;

    galaxyNames = {};
    galaxyData = {};
    index = 1;

    for ii = 1:length(allGalaxyNames)
        for jj = 1:length(myGalaxyNames)
            if strcmp(allGalaxyNames{ii}, myGalaxyNames{jj})
                galaxyNames{index} = allGalaxyNames{ii};
                galaxyData{index} = allGalaxyData{ii};
                index = index + 1;
                break
            end
        end
    end
end

% Get the number of galaxies:
numOfGalaxies = length(galaxyNames);

% Create a cell array for all galaxies:
galaxyMetadata = cell(numOfGalaxies, 1);

% Iterate through all galaxies in order to get rotation curve data:
for ii = 1:numOfGalaxies
    galaxyMetadata{ii}.name = galaxyNames{ii};
    galaxyMetadata{ii}.hubbleType = galaxyData{ii}(1);
end

end


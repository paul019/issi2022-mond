function galaxyNames = filterGalaxies(galaxyNames, galaxyFittingData, chiSquaredReducedThreshold)

galaxyNames_ = {};
index = 1;

for ii = 1:length(galaxyNames)
    if galaxyFittingData{ii}.chiSquaredReducedMin <= chiSquaredReducedThreshold
        galaxyNames_{index} = galaxyNames{ii};
        index = index + 1;
    end
end

galaxyNames = galaxyNames_;

end


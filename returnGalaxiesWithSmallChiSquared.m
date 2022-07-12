function galaxyNames = returnGalaxiesWithSmallChiSquared(galaxies,chiSquaredThreshold)

galaxyNames = {};
index = 1;

for ii = 1:length(galaxies)
    if galaxies{ii}.chiSquaredMin <= chiSquaredThreshold
        galaxyNames{index} = galaxies{ii}.name;
        index = index + 1;
    end
end

end


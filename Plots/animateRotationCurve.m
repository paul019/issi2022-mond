function animateRotationCurve(galaxyName,intFctId,a0,numOfFrames)

fprintf('\n');

path = sprintf('output/anim_%s', galaxyName);

% Create output folder:
if not(isfolder(path))
   mkdir(path)
end

a0Values = [0, logspace(-15,log(a0)/log(10),numOfFrames-1)];

figure
set(gcf, 'Position',  [100, 100, 1000, 750])

for ii = 1:numOfFrames
    mondFit = struct('intFctId',intFctId,'a0',a0Values(ii));

    plotGalaxyRotationCurveWithFits(galaxyName, 0.5, 0.7, {mondFit}, {}, false, false, false);

    title(sprintf('%s, a_0 = %d', galaxyName, a0Values(ii)));

    saveas(gcf,sprintf('%s/frame_%d.png', path, ii));

    hold off

    fprintf('Exported %d of %d frames.\n', ii, numOfFrames);
end

close

fprintf('\n');

end


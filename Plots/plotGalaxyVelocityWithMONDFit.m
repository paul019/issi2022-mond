function plotGalaxyVelocityWithMONDFit(name, a0, MtoLdisk, MtoLbulge, interpolationFunctionIds)
% plot a galaxy named 'name', with Mass to Light ratios of disk MtoLdisk
% (~0.5) and bulge (~0.7).

% Get rotation curve data:
rotationCurveData = ReadRotmodLTGSingle(name);

r = rotationCurveData(:,1);         % in kpc
Vobs = rotationCurveData(:,2);      % in km/s
Vobserr = rotationCurveData(:,3);   % in km/s
Vgas = rotationCurveData(:,4);      % in km/s
Vdisk = rotationCurveData(:,5);     % in km/s
Vbulge = rotationCurveData(:,6);    % in km/s

% Inititlaize variables:
kpcInKm = 3.086*10^16;
velocityLegendArray = {
    'Kinematic data',...
    'Gas velocity',...
    strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')')
};
accelerationLegendArray = {
    'Observed acceleration',...
    'Expected acceleration'
};

% If there are no bulge velocities, set bulgeFlag to false:
if max(Vbulge)==0
    bulgeFlag=false;
else
    bulgeFlag=true;
    velocityLegendArray{end + 1} = strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')');
end

% Add entry to velocityLegendArray:
velocityLegendArray{end + 1} = 'Total baryonic velocity';

% Calculate the baryonic velocity (in km/s):
Vbaryon = sqrt(abs(Vgas).*Vgas+MtoLdisk*abs(Vdisk).*Vdisk+bulgeFlag*MtoLbulge*abs(Vbulge).*Vbulge);

% Calculate corresponding velocities:
Aobs = Vobs.^2 ./ (r*kpcInKm);  % in km/s2
Aobs_min = (Vobs-Vobserr).^2 ./ (r*kpcInKm);  % in km/s2
Aobs_max = (Vobs+Vobserr).^2 ./ (r*kpcInKm);  % in km/s2
Aobserr = (Aobs_max-Aobs_min) ./ 2;   % in km/s2
Aexpected = Vbaryon.^2 ./ (r*kpcInKm);  % in km/s2

% Get the number of interpolation functions:
numOfIntFcts = length(interpolationFunctionIds);

% Initialize cell arrays for MOND fits:
Vmond_r = cell(numOfIntFcts,1);
Vmond   = cell(numOfIntFcts,1);
Amond   = cell(numOfIntFcts,1);

% Calculate the MOND fit for different interpolation functions:
for ii = 1:numOfIntFcts
    intFctId = interpolationFunctionIds{ii};
    intFctName = getInterpolationFunctionName(intFctId);

    [Vmond_r{ii}, Vmond{ii}] = calculateMONDVelocitiesForGalaxy(prepareGalaxyForMOND(name,MtoLdisk,MtoLbulge),a0,intFctId);     % in km/s
    Amond{ii} = Vmond{ii}.^2 ./ (Vmond_r{ii});  % in km/s2

    velocityLegendArray{end + 1} = strcat('MOND fit (', intFctName, ')');
    accelerationLegendArray{end + 1} = strcat('MOND fit (', intFctName, ')');
end

% Prepare figure:
figure

%--------------------------------------------------------------------------
% Subplot 1 (velocities):

subplot(2,1,1);
errorbar(r,Vobs,Vobserr,'.')
hold on;
scatter(r,Vgas)
scatter(r,sqrt(MtoLdisk)*Vdisk)
if bulgeFlag
    scatter(r,sqrt(MtoLbulge)*Vbulge)
end
plot(r,Vbaryon,'--','linewidth',2)
for ii = 1:numOfIntFcts
    plot(Vmond_r{ii}/kpcInKm,Vmond{ii},'--','linewidth',2)
end

title(strcat(name));            
legend(velocityLegendArray, 'Location', 'NorthWest');
grid on;
%set(gca, 'xscale', 'log');
set(gca,'FontSize',15);
%text(1,'FontSize',18);
xlabel 'r [kpc]';
ylabel 'v [km/s]';
axis([0 1.05*max(r) min(Vgas) 1.2*max([max(Vobs),max(Vbaryon)])])

%--------------------------------------------------------------------------
% Subplot 2 (accelerations):

subplot(2,1,2)
hold on;
errorbar(r,Aobs,Aobserr,'.')
plot(r,Aexpected,'--','linewidth',2)
for ii = 1:numOfIntFcts
    plot(Vmond_r{ii}/kpcInKm,Amond{ii},'--','linewidth',2)
end

subplot(2,1,2);
legend(accelerationLegendArray, 'Location', 'NorthWest')
grid on;
set(gca,'FontSize',15);
xlabel 'r [kpc]';
ylabel 'a [km/s^2]';
axis([0 1.05*max(r) min([min(Aobs),min(Aexpected)]) 1.5*max([max(Aobs),max(Aexpected)])])

end

function plotGalaxyVelocityWithFits(name, MtoLdisk, MtoLbulge, mondFits, nfwFits)
% name:         name of the galaxy                  'UGC01281'
% MtoLdisk:     Mass to Light ratio of the disk     0.5
% MtoLbulge:    Mass to Light ratio of the bulge    0.7
% mondFits:     fits using MOND (cell array)        {struct('intFctId','rar','a0',1.2e-13)}
% -> intFctId:  id of the interpolation function
% -> a0:        the parameter a0 (in km/s^2)
% nfwFits:      fits using NFW (cell array)         {struct('p0',3e-23,'R_S',6e16)}
% -> p0:        the parameter p0 (in kg/km^3)
% -> R_S:       the parameter R_S (in km)

% plotGalaxyVelocityWithFits('UGC01281',0.5,0.7,{struct('intFctId','rar','a0',1.2e-13)},{struct('p0',3e-23,'R_S',6e16)});


% Get rotation curve data:
rotationCurveData = ReadRotmodLTGSingle(name);

r = rotationCurveData(:,1);         % in kpc
Vobs = rotationCurveData(:,2);      % in km/s
Vobserr = rotationCurveData(:,3);   % in km/s
Vgas = rotationCurveData(:,4);      % in km/s
Vdisk = rotationCurveData(:,5);     % in km/s
Vbulge = rotationCurveData(:,6);    % in km/s

% Inititlaize variables:
kpcInKm = 3.086e16;
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

%--------------------------------------------------------------------------
% Calculating the MOND fits:

% Get the number of MOND fits:
numOfMONDFits = length(mondFits);

% Initialize cell arrays for MOND fits:
Vmond_r = cell(numOfMONDFits,1);
Vmond   = cell(numOfMONDFits,1);
Amond   = cell(numOfMONDFits,1);

% Prepare rotation curve data for MOND:
rotationCurveData_MOND = prepareGalaxyRotationCurveData(name,MtoLdisk,MtoLbulge,true);

% Calculate the MOND fits:
for ii = 1:numOfMONDFits
    intFctId = mondFits{ii}.intFctId;
    intFctName = getInterpolationFunctionName(intFctId);
    a0 = mondFits{ii}.a0;

    [Vmond_r{ii}, Vmond{ii}] = calculateMONDVelocitiesForGalaxy(rotationCurveData_MOND,a0,intFctId);     % in km/s
    Amond{ii} = Vmond{ii}.^2 ./ (Vmond_r{ii});  % in km/s2

    velocityLegendArray{end + 1} = strcat('MOND fit (', intFctName, ')');
    accelerationLegendArray{end + 1} = strcat('MOND fit (', intFctName, ')');
end

%--------------------------------------------------------------------------
% Calculating the NFW fits:

% Get the number of NFW fits:
numOfNFWFits = length(nfwFits);

% Initialize cell arrays for NFW fits:
Vnfw_r = cell(numOfMONDFits,1);
Vnfw   = cell(numOfMONDFits,1);
Anfw   = cell(numOfMONDFits,1);

% Prepare rotation curve data for NFW:
rotationCurveData_NFW = prepareGalaxyRotationCurveData(name,MtoLdisk,MtoLbulge,false);

% Calculate the NFW fits:
for ii = 1:numOfNFWFits
    p0 = nfwFits{ii}.p0;
    R_S = nfwFits{ii}.R_S;

    [Vnfw_r{ii}, Vnfw{ii}] = calculateNFWVelocitiesForGalaxy(rotationCurveData_NFW,p0,R_S);     % in km/s
    Anfw{ii} = Vnfw{ii}.^2 ./ (Vnfw_r{ii});  % in km/s2

    if numOfNFWFits == 1
        velocityLegendArray{end + 1} = 'NFW fit';
        accelerationLegendArray{end + 1} = 'NFW fit';
    else
        velocityLegendArray{end + 1} = sprintf('NFW fit #%d', ii);
        accelerationLegendArray{end + 1} = sprintf('NFW fit #%d', ii);
    end
end

%--------------------------------------------------------------------------
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
for ii = 1:numOfMONDFits
    plot(Vmond_r{ii}/kpcInKm,Vmond{ii},'--','linewidth',2)
end
for ii = 1:numOfNFWFits
    plot(Vnfw_r{ii}/kpcInKm,Vnfw{ii},'--','linewidth',2)
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
for ii = 1:numOfMONDFits
    plot(Vmond_r{ii}/kpcInKm,Amond{ii},'--','linewidth',2)
end
for ii = 1:numOfNFWFits
    plot(Vnfw_r{ii}/kpcInKm,Anfw{ii},'--','linewidth',2)
end

subplot(2,1,2);
legend(accelerationLegendArray, 'Location', 'NorthWest')
grid on;
set(gca,'FontSize',15);
xlabel 'r [kpc]';
ylabel 'a [km/s^2]';
axis([0 1.05*max(r) min([min(Aobs),min(Aexpected)]) 1.5*max([max(Aobs),max(Aexpected)])])

end

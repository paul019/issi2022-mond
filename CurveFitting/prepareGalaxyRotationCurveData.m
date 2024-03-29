function rotationCurveData = prepareGalaxyRotationCurveData(name,MtoLdisk,MtoLbulge,noImaginaryNumbersFlag)
% prepareGalaxyForMOND

if nargin < 4
    noImaginaryNumbersFlag = true;
end

% Get rotation curve data:
rotationCurveData = ReadRotmodLTGSingle(name);

% rotationCurveData has 8 columns; here is what each column contains:
    % 1 radius in kpc
    % 2 observed velocity in km/s
    % 3 velocity error in km/s
    % 4 gas velocity in km/s
    % 5 disk velocity in km/s
    % 6 bulge velocity in km/s
    % 7 L/pc^2 of disk
    % 8 L/pc^2 of bulge

% Define constants:
kpcInKm = 3.0857e16;

% Define a minimum radius; datapoints with a smaller radius are not
% considered by the code:
minR = 0;  % in kpc

% Adding two columns to rotationCurveData:
    % 9 radius in km
    rotationCurveData(:,9) = rotationCurveData(:,1) * kpcInKm;

    % 10 observed accelreation in km/s^2 (a_obs = v_obs^2 / r)
    rotationCurveData(:,10) = rotationCurveData(:,2).^2./rotationCurveData(:,9);

%{
% Delete rows with data outside the deep-MOND regime:
ii = 1;
dataSize = size(rotationCurveData);
dataLength = dataSize(1);
while ii <= dataLength
    %fprintf('dataLength = %d; ii = %d\n', dataLength, ii);
    insideDeepMONDRegime = rotationCurveData(ii,10) < a0/deepMONDRegimeFactor;

    if (~insideDeepMONDRegime)
        rotationCurveData(ii,:) = [];
        dataLength = dataLength - 1;
    else
        ii = ii + 1;
    end
end
%}

% Delete rows with small radii (see minR):
ii = 1;
dataSize = size(rotationCurveData);
dataLength = dataSize(1);
while ii <= dataLength
    if (rotationCurveData(ii,1) < minR)
        rotationCurveData(ii,:) = [];
        dataLength = dataLength - 1;
    else
        ii = ii + 1;
    end
end

% Calculate the total baryonic velocity:
Vgas = rotationCurveData(:,4);
Vdisk = rotationCurveData(:,5);
Vbulge = rotationCurveData(:,6);

% If there are no bulge velocities, set bulgeFlag to false:
if max(Vbulge)==0
    bulgeFlag=false;
else
    bulgeFlag=true;
end

% Adding a column to rotationCurveData:
    % 11 total baryonic velocity (calculated) in km/s
    rotationCurveData(:,11) = sqrt(abs(Vgas).*Vgas + MtoLdisk*abs(Vdisk).*Vdisk + bulgeFlag*MtoLbulge*abs(Vbulge).*Vbulge);

if noImaginaryNumbersFlag
    % Delete rows where the total baryonic velocity is zero or has an
    % imaginary component:
    ii = 1;
    dataSize = size(rotationCurveData);
    dataLength = dataSize(1);
    while ii <= dataLength
        if (~isreal(rotationCurveData(ii,11))) || rotationCurveData(ii,11) == 0
            rotationCurveData(ii,:) = [];
            dataLength = dataLength - 1;
        else
            ii = ii + 1;
        end
    end
end

end
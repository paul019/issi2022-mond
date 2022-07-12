function [GalsO,GalDataL]=ReadRotmodLTG
% Read Rotmod_LTG files

[GalsO,~]=ReadSPARC;
GalDataL={175};
GalNameL={175};

for ig=1:length(GalsO)
    fname=strcat('Data/Rotmod_LTG/',GalsO{ig},'_rotmod.dat');
    fid = fopen(fname);
    d = fgetl(fid);
    d = fgetl(fid);
    d = fgetl(fid);
    d = fgetl(fid);
    l=1;
    while d~=-1
        b=strfind(d,'	');
        GalNameL{ig}=GalsO{ig}; %reproduce the same name
        GalDataL{ig}(l,1)=str2double(d(1:b(1)-1)); % 1 radius in kpc
        GalDataL{ig}(l,2)=str2double(d(b(1):b(2)-1)); % 2 observed velocity
        GalDataL{ig}(l,3)=str2double(d(b(2):b(3)-1)); % 3 velocity error
        GalDataL{ig}(l,4)=str2double(d(b(3):b(4)-1)); % 4 gas velocity
        GalDataL{ig}(l,5)=str2double(d(b(4):b(5)-1)); % 5 disk velocity
        GalDataL{ig}(l,6)=str2double(d(b(5):b(6)-1)); % 6 bulge velocity
        GalDataL{ig}(l,7)=str2double(d(b(6):b(7)-1)); % 7 L/pc^2 of disk
        GalDataL{ig}(l,8)=str2double(d(b(7):end)); % 8 L/pc^2 of bulge

        l=l+1;
        d = fgetl(fid);
    end
    fclose(fid);
end
GalDataL=transpose(GalDataL);



function [Gals,GalData]=ReadLelliC
% Reads Lelli2016c, a collection of 175 galaxies from the SPARC database.
Ngal=175;

fid = fopen('Data/SPARC_Lelli2016c.txt','r');
for ii=1:98
    d = fgetl(fid);
end
Gals=cell(Ngal,1);
GalData=cell(Ngal,1);
for ii=1:Ngal
    d = fgetl(fid);
    Gals{ii}=strtrim(d(1:12));                         % name
    GalData{ii}(1)=str2num(strtrim(d(13:14)));         % Hubble type
    GalData{ii}(2)=str2num(strtrim(d(17:21)));         % Distance Mpc
    GalData{ii}(3)=str2num(strtrim(d(22:27)));         % Mean error on D Mpc
    GalData{ii}(4)=str2num(strtrim(d(28:30)));         % Distance Method (2)
    GalData{ii}(5)=str2num(strtrim(d(31:35)));         % Inclination deg
    GalData{ii}(6)=str2num(strtrim(d(36:40)));         % Mean error on Inc deg
    GalData{ii}(7)=str2num(strtrim(d(41:48)));         % Total Luminosity at [3.6] 10+9solLum
    GalData{ii}(8)=str2num(strtrim(d(49:56)));         % Mean error on L[3.6] 10+9solLum
    GalData{ii}(9)=str2num(strtrim(d(57:62)));         % Effective Radius at [3.6] kpc
    GalData{ii}(10)=str2num(strtrim(d(63:71)));        % Effective Surface Brightness at [3.6] solLum/pc2
    GalData{ii}(11)=str2num(strtrim(d(72:77)));        % Disk Scale Length at [3.6] kpc, z_d = 0.196*R_d^0.633
    GalData{ii}(12)=str2num(strtrim(d(78:86)));        % Disk Central Surface Brightness at [3.6] solLum/pc2
    GalData{ii}(13)=str2num(strtrim(d(87:94)));        % Total HI mass 10+9 solMass
    GalData{ii}(14)=str2num(strtrim(d(95:100)));       % HI radius at 1 Msun/pc2 - kpc
    GalData{ii}(15)=str2num(strtrim(d(101:106)));      % Asymptotically Flat Rotation Velocity km/s
    GalData{ii}(16)=str2num(strtrim(d(107:112)));      % Mean error on Vflat km/s
    GalData{ii}(17)=str2num(strtrim(d(113:116)));      % Quality Flag (3) 1 = High, 2 = Medium, 3 = Low
end
%% 
% The different ReadXXXXX.m files read the data files of SPARC. Some give
% metadata on the galaxy, some rotation curves. Need to check. In any case,
% it's possible to run plotGalaxyVelocity(...) as below to get a glimpse of
% a rotation curve.
plotGalaxyVelocity('UGC01281',0.5,0.5,true,true) % 0.5 and M/L ratios. true's are flags (see function)
plotGalaxyVelocity('F563-V1',0.5,0.5,true,true) 
%% LelliC: Gals -- galaxy names, GalData -- metadata (see ReadLelliC for information.
% Both structures are array of cells, so to access i'th galaxy, use Gals{i}
% for example
[Gals,GalData]=ReadLelliC;
%% RotmodLTG: GalsO -- galaxy names, GalDataL -- rotation curve data
[GalsO,GalDataL]=ReadRotmodLTG;
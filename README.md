# SPARC analyzer for MOND and NFW

This code was created by the Physics team of the International Summer Science Intsitute (ISSI) 2022 of the Weizmann Institute in Israel. It uses the [open-source SPARC database](http://astroweb.cwru.edu/SPARC/) in order to analyze [galaxy rotation curves with regard to the missing mass problem](https://en.wikipedia.org/wiki/Dark_matter#Galaxy_rotation_curves). The code can perform a [MOND fit](https://en.wikipedia.org/wiki/Modified_Newtonian_dynamics) using six different interpoaltion functions as well as a [Dark Matter NFW fit](https://en.wikipedia.org/wiki/Navarro–Frenk–White_profile).

The following documentation explains the code.

***

## File Structure

- **`Data`**: This folder contains the data from the SPARC database.
- **`BackgroundFunctions`**: This folder contains functions which read data from the SPARC database. Some functions return rotation curves, other functions return galaxy metadata.
- **`OldCode`**: This folder contains code that is not in use anymore.
- **`Plots`**: This folder contains functions which can export the results.
- All functions which are on the main path do the actual computation behind the MOND and NFW fits.

***

## Getting started

- **``**


%% IMPORTANT: Before usage, one has to add the folder "BackgroundFunctions" to the MATLAB path!
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

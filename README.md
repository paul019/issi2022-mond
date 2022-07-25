# SPARC analyser for MOND and NFW

This code was created by the Physics team of the [International Summer Science Intsitute (ISSI)](https://davidson.weizmann.ac.il/en/programs/issi) 2022 of the Weizmann Institute in Israel. It uses the [open-source SPARC database](http://astroweb.cwru.edu/SPARC/) in order to analyse [galaxy rotation curves with regard to the missing mass problem](https://en.wikipedia.org/wiki/Dark_matter#Galaxy_rotation_curves). The code can perform a [MOND fit](https://en.wikipedia.org/wiki/Modified_Newtonian_dynamics) using six different interpoaltion functions as well as a [Dark Matter NFW fit](https://en.wikipedia.org/wiki/Navarro–Frenk–White_profile).

The following documentation explains the code.

***

## File Structure

- **`BackgroundFunctions`**: This folder contains functions which read data from the SPARC database. Some functions return rotation curves, other functions return galaxy metadata.
- **`CurveFitting`**: This folder contains all functions which do the actual computation behind the MOND and NFW fit.
- **`Data`**: This folder contains the data from the SPARC database.
- **`HelpFunctions`**: This folder contains functions which are indirectly used by other functions, e. g. to sort or filter data.
- **`InterpolationFunctions`**: This folder contains functions which regard to the six interpolation functions used for the MOND fit.
- **`OldCode`**: This folder contains code that is not in use anymore.
- **`Output`**: This folder does not exist initially. It is created by the function `exportFits` and used to store all exported results.
- **`Plots`**: This folder contains functions which can export the results.

***

## Installation

1. [Install MATLAB](https://de.mathworks.com/products/matlab.html) on your computer.
2. [Install Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your computer if needed.
3. Open the console on your computer.
    * Use the command `cd` in order to navigate to the directory where you would like this project to be saved.
    * Use the command: `git clone https://github.com/paul019/issi2022-mond`
4. Open MATLAB and navigate to the folder `issi2022-mond` that you just downloaded.
5. Select all subfolders (`BackgroundFunctions`, `CurveFitting`,...) and do: right click &rarr; Add to path &rarr; Selected Folders.
6. Congratulations! The installation is done. You can use the code by running some of the commands below in the MATLAB Command Window.

***

## Getting started

### `exportFits`

#### Input arguments:
| Name | Optional? | Type | Standard value | Meaning |
| ---- | --------- | ---- | -------------- | ------- |
| `galaxyNamesForPlot` | yes | cell array of strings | `{'UGC01281','NGC6503','NGC5055','NGC2366'}` | Galaxy rotation curves for the specified galaxies will be generated. |
| `intFctIds` | yes | cell array of strings | `{ 'linear'; 'standard'; 'exp'; 'rar'; 'simple'; 'toy' }` | The specified interpolation functions will be used for the MOND curve fitting. |
| `minQuality` | yes | integer between 1 and 3 | `3` | Each galaxy in the SPARC database has a quality flag (1 = high, 2 = medium, 3 = low). This input argument determines galaxies of what quality are used for the MOND and NFW fits. |

***

## Credits


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

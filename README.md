# SPARC analyser for MOND and NFW

This code was created by the Physics team of the [International Summer Science Intsitute (ISSI)](https://davidson.weizmann.ac.il/en/programs/issi) 2022 of the Weizmann Institute in Israel. It uses the [open-source SPARC database](http://astroweb.cwru.edu/SPARC/) in order to analyse [galaxy rotation curves with regard to the missing mass problem](https://en.wikipedia.org/wiki/Dark_matter#Galaxy_rotation_curves). The code can perform a [MOND fit](https://en.wikipedia.org/wiki/Modified_Newtonian_dynamics) using six different interpoaltion functions as well as a [Dark Matter NFW fit](https://en.wikipedia.org/wiki/Navarro–Frenk–White_profile).

The following documentation explains the code.

## File Structure

- **`BackgroundFunctions`**: This folder contains functions which read data from the SPARC database. Some functions return rotation curves, other functions return galaxy metadata.
- **`CurveFitting`**: This folder contains all functions which do the actual computation behind the MOND and NFW fit.
- **`Data`**: This folder contains the data from the SPARC database.
- **`HelpFunctions`**: This folder contains functions which are indirectly used by other functions, e. g. to sort or filter data.
- **`InterpolationFunctions`**: This folder contains functions which regard to the six interpolation functions used for the MOND fit.
- **`OldCode`**: This folder contains code that is not in use anymore.
- **`Output`**: This folder does not exist initially. It is created by the function `exportFits` and used to store all exported results.
- **`Plots`**: This folder contains functions which can export the results.

## Installation

1. [Install MATLAB](https://de.mathworks.com/products/matlab.html) on your computer.
2. [Install Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your computer if needed.
3. Open the console on your computer.
    * Use the command `cd` in order to navigate to the directory where you would like this project to be saved.
    * Use the command: `git clone https://github.com/paul019/issi2022-mond`
4. Open MATLAB and navigate to the folder `issi2022-mond` that you just downloaded.
5. Select all subfolders (`BackgroundFunctions`, `CurveFitting`,...) and do: right click &rarr; Add to path &rarr; Selected Folders.
6. Congratulations! The installation is done. You can use the code by running some of the commands below in the MATLAB Command Window.

## Getting started

### `exportFits`

#### What this function does:

This function is built on top of nearly every other function. It performs the specified MOND and NFW fits and exports the obtained data in the form of several tables and plots. Note that this function can take a while to run as it has to do some heavy computing.

#### Input arguments:

| Name | Optional? | Type | Standard value | Meaning |
| ---- | --------- | ---- | -------------- | ------- |
| `galaxyNamesForPlot` | yes | cell array of strings | `{'UGC01281','NGC6503','NGC5055','NGC2366'}` | Galaxy rotation curves for the specified galaxies will be generated. However, more galaxies will be used to do the MOND and NFW fits (see `minQuality`). |
| `intFctIds` | yes | cell array of strings | `{ 'linear'; 'standard'; 'exp'; 'rar'; 'simple'; 'toy' }` | The specified interpolation functions will be used for the MOND curve fitting. |
| `minQuality` | yes | integer between 1 and 3 | `3` | Each galaxy in the SPARC database has a quality flag (1 = high, 2 = medium, 3 = low). All galaxies from the SPARC database which comply with the minimum quality will be used to do the MOND and NFW fits. |

#### Output arguments:

| Name | Type | Meaning |
| ---- | ---- | ------- |
| `fitIds` | cell array of strings | In this array, there is one cell for every fit that was done, i. e. each MOND interpolation function plus one cell for the NFW fit. |
| `galaxyFittingDataArray` | cell array of cell arrays of structs | In this array, there is one cell for every fit (see `fitIds`). Each cell contains another cell array; here, there is one cell for every galaxy containing details of the fit. |

Additional to these outputs, the function also writes `png` and `csv` files to the `Output` folder:

| File name | Content |
| --------- | ------- |
| `fits_overview.csv` | A table with general information on the MOND and NFW fits. |
| `fits_galaxy_details.csv` | A table with detailed information on the MOND and NFW fits for every galaxy. |
| `intFcts.png` | A plot of all used interpolation functions. |
| `MONDfits_MSWD_vs_a0.png` | A plot of the [Mean Square Weighted Deviation (MSWD)](https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic) per degree of freedom against different values of a<sub>0</sub> for each used MOND interpolation function. |
| `fits_histogram.png` | A figure with one plot for each fitting method, i. e. each MOND interpolation function and NFW. Each plot shows a histogramm of the MSWD per degree of freedom for the respective fit. |
| `GALAXYNAME_v_vs_r.png` | A plot showing the rotation curves for the respective galaxy. The rotation curves include observed velocities, expected Newtonian velocities, and the fitted curves for MOND and NFW. Note that the MOND fits rely on the best a<sub>0</sub> value overall of the respective interpolation function. |

#### Example usage:

`exportFits;`

 ***

### `animateRotationCurve`

#### What this function does:

This function returns an animation that shows how the expected Newtonian velocity curve transitions into the MOND curve smoothly.

#### Input arguments:

| Name | Optional? | Type | Meaning |
| ---- | --------- | ---- | ------- |
| `galaxyName` | no | string | The galaxy for which the animation is created. |
| `intFctId` | no | string | The interpoaltion function which is used for the animation. |
| `a0` | no | floating point number | The a<sub>0</sub> value which is used for the animation (in km/s<sup>2</sup>). |
| `numOfFrames` | no | integer | The number of animation frames. |

#### Output arguments:

None. The function creates a folder with the path `Output/anim_GALAXYNAME` and saves one `png` file for every animation frame in there.

#### Example usage:

`animateRotationCurve('NGC6503','rar',1.2e-13,20);`

## Credits

This project was done by Andrés, Christian, Shashank, Nimrod, Paul, Aamod, and Yin Lam. Thanks to Abishek for mentoring us throughout the program!

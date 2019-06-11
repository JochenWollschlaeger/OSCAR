The MATLAB (The MathWorks, USA) code presented here is for calculating absorption coefficient spectra based on the light intensity spectra provided by the Online Hyperspectral Integrating Cavity Absorption Meter (OSCAR, TriOS GmbH, Germany).

Please cite the software using the following doi:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3243405.svg)](https://doi.org/10.5281/zenodo.3243405)

Instructions for software use:

1.Create a folder "Rawdata" in the same folder where the scripts are.

2.Put the "RAW_DARK.dat" and the corresponding "RAW_LIGH.dat" files, which are in the .tar file downloaded from the OSCAR in that folder. The files can also be renamed.

3.Fill out the "User_Input.xlsx" spreadsheet and put it in the same location where also the scripts are.

4.Read the instructions in script "A_Import_OSCAR_Data", make the appropriate changes in the "Parameters" section, and execute the script.

5.Check the result files and make changes, if required (e.g. delete wrong measurements or give temperature and salinity information in the appropriate columns).

6.Read the instructions in script "B_Calculate_Reflectivity", make the appropriate changes in the "Parameters" section, and execute the script.

7.Check the result file and make changes, if necessary.

8.Read the instructions in script "C_Calculate_Absorption", make the appropriate changes in the "Parameters" section, and execute the script.

The MATLAB (The MathWorks, USA) code presented here is for calculating absorption coefficient spectra based on the light intensity spectra provided by the Online Hyperspectral Integrating Cavity Absorption Meter (OSCAR, TriOS GmbH, Germany).

Short instructions for software use:

1. Unzip "OSCAR_Code.zip" to retrieve the MATLAB scripts, the "Subfunctions" folder, and the "User_Input_Template.xlsx" file.

2. Create a folder "Rawdata" in the same folder where the scripts are.

3. Put the "RAW_DARK" and the corresponding "RAW_LIGH" files extracted from the .tar file downloaded from the OSCAR in that folder (either as .dat or .csv). Renaming of files is possible.

4. Rename the "User_Input_Template.xlsx" into "User_Input.xlsx".

5. Fill out the spreadsheet with the necessary information according to the instructions given in the file and put it in the same location where the scripts are.

6. Read the instructions in script "A_Import_OSCAR_Data", make the appropriate changes in the "Parameters" section, and execute the script.

7. Check the result files and make changes, if required (e.g. delete wrong measurements or give temperature and salinity information in the appropriate columns).

8. Read the instructions in script "B_Calculate_Reflectivity", make the appropriate changes in the "Parameters" section, and execute the script.

9. Check the result file and make changes, if necessary.

10. Read the instructions in script "C_Calculate_Absorption", make the appropriate changes in the "Parameters" section, and execute the script.

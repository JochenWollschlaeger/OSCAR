% Script for calculating the reflectivity of an OSCAR integrating cavity
% (TriOS GmbH, Germany) 
%
% The script requires the output of the "A_Import_OSCAR.m" script,
% specifically the "CalibrationData.txt" and the "Constants.txt".
%
% The result of the script is a textfile ("Reflectivity.txt") that contains
% the reflectivity spectra as well as the begin and end of each calibration.
%
% The reflectivity spectra should be checked for their reasonability, and
% incorrect calibrations can be deleted before calculating the absorption
% coefficients of the sample measurements. Do not change timestamp formats
% or length of individual vectors. This will cause problems for subsequent
% scripts.
%
% Jochen Wollschläger
% Institute for Chemistry and Biology of the Marine Environment
% University of Oldenburg
% Contact: jochen.wollschlaeger@uol.de
% 
% Changelog:
%
% Version 1.0: Finished (07.06.2019)
%
% License information:
% Copyright (c) [2019] [Jochen Wollschläger]
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
% OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
% THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%--------------------------------------------------------------------------

clc
clear
BasicPath=pwd;

%% Parameters

%Path of the folder that contains the output of the previously executed 
%script
InputDataFolder=[BasicPath,'\Results'];

%Name of the folder that will contain the output files
%If not existing, the folder will be created
ResultsFolder=[BasicPath,'\Results'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FROM HERE, ONLY MAKE MODIFICATIONS IF YOU KNOW WHAT YOU ARE DOING!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if all folders and files are available

cd(BasicPath)

%Output from previous script(s)
if exist(InputDataFolder,'dir')==7
    cd(InputDataFolder)
    if ~exist('CalibrationData.txt','file')==2
        error('"CalibrationData.txt" is missing!')
    end
    if ~exist('Constants.txt','file')==2
        error('"Constants.txt" is missing!')
    end
else
    error(['The folder ',InputDataFolder,' is missing!'])
end

cd(BasicPath)

%% Data import

cd(InputDataFolder)

%Calibration data
fid=fopen('CalibrationData.txt','r');
Data=textscan(fid,[repmat('%s', 1, 3000) '%*[^\n]'], 'delimiter', '\t', 'collectoutput', true);
Data=[Data{:}];
fclose(fid);

DesiredFormat=Data{1,3};
DesiredFormat=strsplit(DesiredFormat,{'[',']'});
DesiredFormat=DesiredFormat{2};

Data=Data(2:end,:);

Calib_Name=unique(Data(:,1),'stable');

AbsCoeff_Standard=Data(ismember(Data(:,2),{'AbsCoeff_Standard'}),7:end);
AbsCoeff_Standard=str2double(AbsCoeff_Standard);
AbsCoeff_Standard=AbsCoeff_Standard(:,~isnan(AbsCoeff_Standard(1,:)));

Calib_Ref_Timestamp_Start=Data(ismember(Data(:,2),{'Calib_Ref'}),3);
Calib_Ref_Temp=Data(ismember(Data(:,2),{'Calib_Ref'}),5);
Calib_Ref_Temp=str2double(Calib_Ref_Temp);
Calib_Ref_Data=Data(ismember(Data(:,2),{'Calib_Ref'}),7:end);
Calib_Ref_Data=str2double(Calib_Ref_Data);
Calib_Ref_Data=Calib_Ref_Data(:,~isnan(Calib_Ref_Data(1,:)));

Calib_Standard_Timestamp_End=Data(ismember(Data(:,2),{'Calib_Standard'}),4);
Calib_Standard_Temp=Data(ismember(Data(:,2),{'Calib_Standard'}),5);
Calib_Standard_Temp=str2double(Calib_Standard_Temp);
Calib_Standard_Data=Data(ismember(Data(:,2),{'Calib_Standard'}),7:end);
Calib_Standard_Data=str2double(Calib_Standard_Data);
Calib_Standard_Data=Calib_Standard_Data(:,~isnan(Calib_Standard_Data(1,:)));
clear Data ans fid

%Constants data
fid=fopen('Constants.txt','r');
Data1=textscan(fid,'%s',3,'Delimiter','\t');
Data1=[Data1{:}];
Data2=textscan(fid,[repmat('%s', 1, 3000) '%*[^\n]'], 'delimiter', '\t', 'collectoutput', true);
Data2=[Data2{:}];
fclose(fid);

Wavelengths=str2double(Data2(1,2:end));
AbsCoeff_Water=str2double(Data2(2,2:end));
TempCoeff=str2double(Data2(3,2:end));

Wavelengths=Wavelengths(~isnan(Wavelengths));
AbsCoeff_Water=AbsCoeff_Water(~isnan(AbsCoeff_Water));
TempCoeff=TempCoeff(~isnan(TempCoeff));

r=strsplit(Data1{1,1});
r=str2double(r(end));
r0=strsplit(Data1{2,1});
r0=str2double(r0(end));
AbsCoeff_WaterT=strsplit(Data1{3,1});
AbsCoeff_WaterT=str2double(AbsCoeff_WaterT(end));
clear Data1 Data2 fid ans

cd(BasicPath)

%% Creating NaN vectors in case no calibration measurements have been performed

if isempty(AbsCoeff_Standard)
    AbsCoeff_Standard=NaN(size(Wavelengths));
    Calib_Ref_Data=NaN(size(Wavelengths));
    Calib_Standard_Data=NaN(size(Wavelengths));
end

%% Check for available temperature

TempCorr=cell(size(Calib_Name));
for i=1:length(Calib_Name)
    if isnan(Calib_Ref_Temp(i)) || isnan(Calib_Standard_Temp(i))
        TempCorr(i)={'No'};
        Calib_Ref_Temp(i)=AbsCoeff_WaterT;
        Calib_Standard_Temp(i)=AbsCoeff_WaterT;
    else
        TempCorr(i)={'Yes'};
    end
end
clear i

%% Calculation of reflectivity

Reflectivity=zeros(size(AbsCoeff_Standard));
for i=1:size(Reflectivity,1)
    Reflectivity(i,:)=reflectivity_calc(Calib_Ref_Data(i,:),...
                                        Calib_Standard_Data(i,:),...
                                        Calib_Ref_Temp(i),...
                                        Calib_Standard_Temp(i),...
                                        AbsCoeff_Standard(i,:),...
                                        AbsCoeff_Water,...
                                        AbsCoeff_WaterT,...
                                        TempCoeff,...
                                        r,r0);
end
clear i

%% Saving the data

%Going to or creating the folder where the working data will be saved
if exist(ResultsFolder,'dir')==0
    mkdir(ResultsFolder)
end

cd(ResultsFolder)

SaveFileName='Reflectivity.txt';

HeaderFormat=['%s\t%s\t%s\t%s',repmat('\t%d',1,size(Wavelengths,2)),'\r\n'];
Header=[{'Calibration No./Name'} {['Timestamp Start [',DesiredFormat,']']} {['Timestamp End [',DesiredFormat,']']} {'Temperature correction'} num2cell(Wavelengths)];
DataFormat=['%s\t%s\t%s\t%s',repmat('\t%f',1,size(Wavelengths,2)),'\r\n'];
Data=[Calib_Name Calib_Ref_Timestamp_Start Calib_Standard_Timestamp_End TempCorr num2cell(Reflectivity)];

fid=fopen(SaveFileName,'w');
fprintf(fid,HeaderFormat,Header{1,:});
for i=1:size(Data,1)
    fprintf(fid,DataFormat,Data{i,:});
end
fclose(fid);
clear i fid ans SaveFileName Header HeaderFormat Data DataFormat

cd(BasicPath)

Properties=struct('WindowStyle','non-modal','Interpreter','tex');
msgbox('\bf \fontsize{12} \fontname{Arial} Script finished',Properties)
clear Properties

%% Non build-in functions required from the script

%reflectivity_calc
function [rho]=reflectivity_calc(ref_data,std_data,ref_temp,std_temp,abs_std,abs_water,abs_waterT,tempcoeff,r,r0)
% This function calculates the absorption coefficients of the optically 
% active constituents in a sample inside an integrating cavity.
% The equations are adopted from Röttgers et al (2005) Practical test of a 
% point-source integrating cavity absorption meter: the performance of 
% different collector assemblies. Appl. Opt. 44, 5549–5560.
%
% Input parameters:
% ref_data (vector): Reference spectrum
% std_data (vector): Standard spectrum
% ref_temp (scalar): Temperature of the reference measurement
% std_temp (scalar): Temperature of the standard measurement
% abs_std (vector): Absorption coefficient spectrum of the standard used
%                   during the calibration
% abs_water (vector): Absorption coefficient spectrum of pure water
% temp_coeff (vector): Coefficients for correction of the absorption 
%                       coefficients of pure water for temperature effects
% r (scalar): Radius of the cavity
% r0 (scalar): Radius of the cavity minus radius of the point light source
%
% Output parameter (data format):
% rho (vector): Reflectivity spectrum inside the cavity
%
% Jochen Wollschläger
% Institute for Chemistry and Biology of the Marine Environment
% University of Oldenburg
% Contact: jochen.wollschlaeger@uol.de
% 
% Changelog:
%
% Version 1.0: Finished 18.01.2019
%--------------------------------------------------------------------------

%Calculating the absorption of the used solutions 
%(including temperature correction)
Abs_Ref=abs_water+((ref_temp-abs_waterT).*tempcoeff);
Abs_Std=abs_std+abs_water+((std_temp-abs_waterT).*tempcoeff);

%Calculating_transmission
Trans=std_data./ref_data;
 
%Calculating reflectivity (rho)
PsaAr=(1-exp(-2*Abs_Std*r).*(2*Abs_Std*r+1))./(2*Abs_Std.^2*r^2);
PsaBr=(1-exp(-2*Abs_Ref*r).*(2*Abs_Ref*r+1))./(2*Abs_Ref.^2*r^2);
aAr0=Abs_Std*r0;  
aBr0=Abs_Ref*r0;

rho=(Trans.*exp(-aBr0).*PsaBr-exp(-aAr0).*PsaAr)./...
    ((Trans.*exp(-aBr0).*PsaAr.*PsaBr)-(exp(-aAr0).*PsaBr.*PsaAr));

end


% Script for calculating absorption coefficient spectra from measurement
% with an Online Hyperspectral Integrating Cavity Absorption Meter OSCAR
% (TriOS GmbH, Germany)
%
% The script requires the output of the scripts "A_Import_OSCAR.m" and 
% "B_Calculate_Reflectivity.m", specifically "Constants.txt",
% "ReferenceData.txt", "SampleData.txt", and "Reflectivity.txt".
%
% The result of this script is a textfile ("Absorption.txt"). It contains
% the absorption coefficient spectra based calculated from the difference 
% of the light intensity during the reference and sample measurements,
% taking into account the reflectivity of the cavity which determines the 
% optical path length. It gives also a quality flag that provides a first
% guess whether the spectrum is of reasonable quality or not (1: quality 
% ok, 0: quality questionable).
%
% For the time periods in between the individual calibration and reference
% measurements, respectively, the data are interpolated to provide
% infomation for the absorption coefficient calculation. The interpolation
% method can be defined in the "Parameters" section below (see also the
% documention regarding the interp1.m function).
% 
% Besides a correction for temperature, the absorption coefficient values 
% also have to be corrected for the influence of salinity on the pure water
% absorption. In case these data are not available, a default salinity
% value can be specified in the "Parameters" section that is used unless no
% other salinity values are given in the "SampleData.txt" file.
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
%scripts
InputDataFolder=[BasicPath,'\Results'];

%Name of the folder that will contain the output files
%If not existing, the folder will be created
ResultsFolder=[BasicPath,'\Results'];

%Salinity value that is used for salinity correction unless no other value
%is specified
DefaultSalinity=0;

%Interpolation method for the reference and reflectivity data
InterpMethod='linear'; % other options are 'previous' or 'spline'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FROM HERE, ONLY MAKE MODIFICATIONS IF YOU KNOW WHAT YOU ARE DOING!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if all folders and files are available

cd(BasicPath)

%Output from previous script(s)
if exist(InputDataFolder,'dir')==7
    cd(InputDataFolder)
    if ~exist('Constants.txt','file')==2
        error('"Constants.txt" is missing!')
    end
    if ~exist('ReferenceData.txt','file')==2
        error('"ReferenceData.txt" is missing!')
    end
    if ~exist('SampleData.txt','file')==2
        error('"SampleData.txt" is missing!')
    end
    if ~exist('Reflectivity.txt','file')==2
        error('"Reflectivity.txt" is missing!')
    end
else
    error(['The folder ',FolderInOutName,' is missing!'])
end

%% Import reflectivity data

cd(InputDataFolder)

fid=fopen('Reflectivity.txt','r');
Data=textscan(fid,[repmat('%s', 1, 3000) '%*[^\n]'], 'delimiter', '\t', 'collectoutput', true);
Data=[Data{:}];
fclose(fid);

DesiredFormat=Data{1,2};
DesiredFormat=strsplit(DesiredFormat,{'[',']'});
DesiredFormat=DesiredFormat{2};

Data=Data(2:end,:);

DesiredFormat_Used=DesiredFormat(1:length(Data{1,2}));

Reflectivity=Data(:,5:end);
Reflectivity=Reflectivity(:,~cellfun(@isempty,Reflectivity(1,:)));
Reflectivity=str2double(Reflectivity);

Reflectivity_Timestamp_Start=datenum(Data(:,2),DesiredFormat_Used);
Reflectivity_Timestamp_End=datenum(Data(:,3),DesiredFormat_Used);

Reflectivity_Tempcorr=Data(:,4);
clear Data fid ans DesiredFormat_Used

%% Import reference data

fid=fopen('ReferenceData.txt','r');
Data=textscan(fid,[repmat('%s', 1, 3000) '%*[^\n]'], 'delimiter', '\t', 'collectoutput', true);
Data=[Data{:}];
fclose(fid);

Data=Data(2:end,:);

DesiredFormat_Used=DesiredFormat(1:length(Data{1,2}));

Ref_Timestamp_Start=datenum(Data(:,2),DesiredFormat_Used);
Ref_Timestamp_End=datenum(Data(:,3),DesiredFormat_Used);
Ref_Temp=str2double(Data(:,4));
Ref_Data=Data(:,6:end);
Ref_Data=Ref_Data(:,~cellfun(@isempty,Ref_Data(1,:)));
Ref_Data=str2double(Ref_Data);
clear Data fid ans DesiredFormat_Used

%% Import sample data

%Obtain the file format from the file
%(increases the import speed of large data amounts)
fid=fopen('SampleData.txt','r');
Data=textscan(fid,[repmat('%s', 1, 3000) '%*[^\n]'],1, 'delimiter', '\t', 'collectoutput', true);
Data=[Data{:}];
Data=Data(:,~cellfun(@isempty,Data(1,:)));
IndexOfNumerics=~isnan(str2double(Data));
FileFormat=[repmat('%s',1,sum(~IndexOfNumerics)),repmat('%f',1,sum(IndexOfNumerics))];
fclose(fid);
clear IndexOfNumerics Data fid ans

%Import the data with the obtained file format
fid=fopen('SampleData.txt','r');
Data=textscan(fid,[FileFormat '%*[^\n]'], 'delimiter', '\t', 'collectoutput', true);
fclose(fid);
Sample_Name=Data{1,1}(2:end,1);
Sample_Timestamp=datenum(Data{1,1}(2:end,2),DesiredFormat);
Sample_Temp=str2double(Data{1,1}(2:end,3));
Sample_Salt=str2double(Data{1,1}(2:end,4));
Sample_Data=Data{1,2}(2:end,:);
clear Data fid ans FileFormat

%% Import Constants data

fid=fopen('Constants.txt','r');
Data1=textscan(fid,'%s',3,'Delimiter','\t');
Data1=[Data1{:}];
Data2=textscan(fid,[repmat('%s', 1, 3000) '%*[^\n]'], 'delimiter', '\t', 'collectoutput', true);
Data2=[Data2{:}];
fclose(fid);

Wavelengths=str2double(Data2(1,2:end));
AbsCoeff_Water=str2double(Data2(2,2:end));
TempCoeff=str2double(Data2(3,2:end));
SaltCoeff=str2double(Data2(4,2:end));

Wavelengths=Wavelengths(~isnan(Wavelengths));
AbsCoeff_Water=AbsCoeff_Water(~isnan(AbsCoeff_Water));
TempCoeff=TempCoeff(~isnan(TempCoeff));
SaltCoeff=SaltCoeff(~isnan(SaltCoeff));

r=strsplit(Data1{1,1});
r=str2double(r(end));
r0=strsplit(Data1{2,1});
r0=str2double(r0(end));
AbsCoeff_WaterT=strsplit(Data1{3,1});
AbsCoeff_WaterT=str2double(AbsCoeff_WaterT(end));
clear Data1 Data2 fid ans

cd(BasicPath)

%% Assigning reflectivity and reference data to the sample measurements

%Interpolation of data over all timestamps
AllTimestamps=unique([Reflectivity_Timestamp_Start;Reflectivity_Timestamp_End;Ref_Timestamp_Start;Ref_Timestamp_End;Sample_Timestamp]);
Ref_Data_interp=interp1([Ref_Timestamp_Start;Ref_Timestamp_End],[Ref_Data;Ref_Data],AllTimestamps,InterpMethod,'extrap');
Ref_Temp_interp=interp1([Ref_Timestamp_Start;Ref_Timestamp_End],[Ref_Temp;Ref_Temp],AllTimestamps,InterpMethod,'extrap');
Reflectivity_interp=interp1([Reflectivity_Timestamp_Start;Reflectivity_Timestamp_End],[Reflectivity;Reflectivity],AllTimestamps,InterpMethod,'extrap');

Reflectivity_Tempcorr_interp=cell(size(AllTimestamps));
Timestamps_Working=[Reflectivity_Timestamp_End;AllTimestamps(end)];
for i=1:size(Timestamps_Working,1)-1
    Index=[AllTimestamps>=Timestamps_Working(i) AllTimestamps<=Timestamps_Working(i+1)];
    Index=sum(Index,2)>0;
    Reflectivity_Tempcorr_interp(Index)=Reflectivity_Tempcorr(i);
end
clear Index i Timestamps_Working

%Selecting the data for the timestamps of the sample measurements
Index=ismember(AllTimestamps,Sample_Timestamp);
Ref_Data=Ref_Data_interp(Index,:);
Ref_Temp=Ref_Temp_interp(Index,:);
Reflectivity=Reflectivity_interp(Index,:);
Reflectivity_Tempcorr=Reflectivity_Tempcorr_interp(Index,:);
clear Index Ref_Mean_interp Ref_Temp_interp Reflectivity_interp...
      Reflectivity_Tempcorr_interp

%% Check for available temperature and salinity data

TempCorr=cell(size(Sample_Data,1),1);
for i=1:size(Sample_Data,1)
    if isnan(Sample_Temp(i))||isnan(Ref_Temp(i))
        TempCorr(i)={'No'};
        Sample_Temp(i)=20.1;
        Ref_Temp(i)=20.1;
    else
        TempCorr(i)={'Yes'};
    end
end
clear i

SaltCorr=cell(size(Sample_Data,1),1);
for i=1:size(Sample_Data,1)
    if isnan(Sample_Salt(i))
        SaltCorr(i)={['Used default salinity (',num2str(DefaultSalinity),')']};
        Sample_Salt(i)=DefaultSalinity;
    else
        SaltCorr(i)={'Yes'};
    end
end
clear i
clear DefaultSalinity

%% Absorption coefficient calculation
AbsorptionCoeff=zeros(size(Sample_Data));
QualityFlag=NaN(size(Sample_Data,1),1);
AbsorptionSpectraToCalculate=size(AbsorptionCoeff,1);
for i=1:size(AbsorptionCoeff,1)
    [abs_coeff,quality]=absorption_calc(Ref_Data(i,:),...
                                        Sample_Data(i,:),...
                                        Reflectivity(i,:),...
                                        AbsCoeff_Water,...
                                        AbsCoeff_WaterT,...
                                        Ref_Temp(i),...
                                        Sample_Temp(i),...
                                        TempCoeff,...
                                        Sample_Salt(i),...
                                        SaltCoeff,...
                                        r,r0,...
                                        Wavelengths);
    AbsorptionCoeff(i,:)=abs_coeff;
    QualityFlag(i)=quality;
    AbsorptionSpectraToCalculate=AbsorptionSpectraToCalculate-1;
    display(AbsorptionSpectraToCalculate)
end
clear i abs_coeff quality AbsorptionSpectraToCalculate Ref_Mean_Used...
      Ref_Temp_Used Reflectivity_Used Sample Sample_Salt Sample_Temp

%% Convert the timestamp in readable format

Sample_Timestamp=cellstr(datestr(Sample_Timestamp,DesiredFormat));

%% Save the data

%Going to or creating the folder where the working data will be saved
if exist(ResultsFolder,'dir')==0
    mkdir(ResultsFolder)
end

cd(ResultsFolder)

SaveFileName='Absorption.txt';

HeaderFormat=['%s\t%s\t%s\t%s\t%s\t%s',repmat('\t%d',1,size(Wavelengths,2)),'\r\n'];
Header=[{'Sample No./Name'} {['Timestamp [',DesiredFormat,']']} {'Reflec. Temp. corr.'} {'Temp. corr.'} {'Sal. corr.'} {'QualityFlag'} num2cell(Wavelengths)];
DataFormat=['%s\t%s\t%s\t%s\t%s\t%.0f',repmat('\t%f',1,size(Wavelengths,2)),'\r\n'];
Data=[Sample_Name Sample_Timestamp Reflectivity_Tempcorr TempCorr SaltCorr num2cell(QualityFlag) num2cell(AbsorptionCoeff)];

fid=fopen(SaveFileName,'w');
fprintf(fid,'%s\r\n',['Interpolation method for reference and reflectivity: ',InterpMethod]);
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

%absorption_calc
function [abs_coeff,quality]=absorption_calc(ref_data,sample_data,rho,abs_water,abs_waterT,ref_temp,sample_temp,tempcoeff,salinity,saltcoeff,r,r0,wavelengths)
% This function calculates the absorption coefficients of the optically 
% active constituents in a sample inside an integrating cavity.
% The equations are adopted from Röttgers et al (2005) Practical test of a 
% point-source integrating cavity absorption meter: the performance of 
% different collector assemblies. Appl. Opt. 44, 5549–5560.
%
% Input:
% ref_data (vector): Reference spectrum
% sample_data (vector): Sample spectrum
% rho (vector): Reflectivity spectrum inside the cavity
% abs_water (vector): Absorption coefficient spectrum of pure water
% ref_temp (scalar): Temperature of the reference measurement
% sample_temp (scalar): Temperature of the sample measurement
% tempcoeff (vector): Coefficients for correction of the pure water 
%                     absorption for the effects of temperature
% salinity (scalar): Salinity of the sample
% saltcoeff (vector): Coefficients for correction of the pure water 
%                     absorption for the effects of salinity
% r (scalar): Radius of the cavity
% r0 (scalar): Radius of the cavity minus radius of the point light source
% wavelengths (vector): Wavelength range covered
%
% Output:
% abs_coeff (vector): Wavelength-specific absorption coefficients
% quality (scalar): First guess quality flag (1 = probably ok, 
%                   0 = potentially bad) based on the following criteria: 
%                   The absorption coefficient spectrum has to be positive
%                   at 700 nm and the slope between 700 and 704 nm has to
%                   be smaller than zero.
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

%Calculation of transmission
trans=sample_data./ref_data;
trans=loess_m(1:length(trans),trans,1:length(trans),0.05);

%Calculating water absorption of the reference (correction for temperature)
abs_wass_ref=abs_water+((ref_temp-abs_waterT)*tempcoeff);

%Calculating water absorption of the sample (correction for temperature and
%salinity)
abs_wass_sample=abs_water+((sample_temp-abs_waterT)*tempcoeff+salinity*saltcoeff); 

%Calculating absorption of sample
abr=r*abs_wass_ref; 
Psab=(1-exp(-2*abr).*(2*abr+1))./(2*abr.^2);
abs_coeff=zeros(size(trans));
for i=1:size(trans,2) %Calculation of absorption per wavelength
        [x,~,~,~]=fminsearch(@absfunc,abs_wass_sample(i),[],trans(i),...
                             rho(i),Psab(i),r0,abs_wass_ref(i),r);
        abs_coeff(i)=x-abs_wass_sample(i);
end

%Smoothing of the resulting absorption coefficient spectrum
abs_coeff=loess_m(1:length(abs_coeff),abs_coeff,1:length(abs_coeff),0.05);

%Quality check
if max(wavelengths)<700
    quality=2;
else    
    Slope=(abs_coeff(wavelengths==650)-abs_coeff(wavelengths==690))/(650-690);
    if abs_coeff(wavelengths==700)>0 && Slope<0
        quality=1;
    else
        quality=3;
    end
end

end

%absfunc
function diff=absfunc(a,T,rho,Psaref,r0,aref,r)
% This function is used in the iterational determination of the absorption
% coefficient of the absorption_calc.m function.
% It calculates the difference between the measured and estimated
% transmission. The latter is a function of the absorption.
% The calculations are based on equations given in 
% Röttgers et al (2005) Practical test of a point-source integrating cavity
% absorption meter: the performance of different collector assemblies. 
% Appl. Opt. 44, 5549–5560.
%
% Input:
% a (scalar): Absorption coefficient at a specific wavelength
% T (scalar): Transmission at a specific wavelength
% rho (scalar): Reflectivity of the cavity at a specific wavelength
% Psaref (scalar): Probability of a photon reaching the cavity wall
% r0 (scalar): Radius of the cavity minus radius of the point light source
% aref (scalar): Absorption coefficient at a specific wavelength
% r (scalar): Radius of the cavity
%
% Output:
% diff (scalar): Squared difference between measured and estimated
%                transmission.
%
% Original code provided by Rüdiger Röttgers
% Institute for Coastal Research
% Helmholtz Zentrum Geesthacht.
%
% Modified by Jochen Wollschläger
% Institute for Chemistry and Biology of the Marine Environment
% University of Oldenburg
% Contact: jochen.wollschlaeger@uol.de
%
% Changelog:
%
% Version 1.0: Finished 18.01.2019
%--------------------------------------------------------------------------

Psa=(1-exp(-2*a*r)*(2*a*r+1))/(2*a^2*r^2);

tguess=exp(-r0*(a-aref))*(1-rho*Psaref)/(1-rho*Psa)*Psa/Psaref;

diff=(T-tguess)*(T-tguess);
end

%loess_m
function yi=loess_m(xr,yr,xi,alpha)
% loess (low S)-type of interpolation
% B. Efron; An Introduction to the bootstrap; ISBN 0-412-04231-2; p.77
% 
% Input:
% xr (vector): x values of the original vector
% yr (vector): y values of the original vector
% xi (vector): New x values, where the interpolation is to be calculated
% alpha (scalar): Relative amount of points to be used in local smoothing
%                 (ranging from 0-1)
% 
% Output:
% yi (vector): Interpolated data
% 
% Original code provided by Rüdiger Röttgers
% Institute for Coastal Research
% Helmholtz Zentrum Geesthacht.
%
% Modified by Jochen Wollschläger
% Institute for Chemistry and Biology of the Marine Environment
% University of Oldenburg
% Contact: jochen.wollschlaeger@uol.de
%
% yr can now also be a matrix (data organized in rows). Allows multiple 
% interpolations at once.
% 
% Changelog:
%
% Version 1.0: Finished 18.01.2019
%--------------------------------------------------------------------------

if size(xr,2) ~= size(yr,2) 
    fprintf('loess_m expects equal lengths for first two input variables\n');
    return
end
if (alpha <= 0.) || (alpha >= 1.) 
    fprintf('loess_m expects alpha in (0,1) but finds alpha=%f\n',alpha);
    return
end

np=round(alpha*length(xr));
[~,h]=sort(xr);
x=xr(h);
y=yr(:,h);
yi=zeros(size(yr,1),length(xi));
for i=1:length(xi)
    [~,h]=sort(abs(x-xi(i)));
    h=h(1:np);
    xt=x(h);
    yt=y(:,h);
    u=abs(xt-xi(i));
    u=u/max(u);
    w=(1.-u.^3).^3;
    [a,b]=fitw_m(xt,yt,w);
    yi(:,i)=a+b.*xi(i);
end
end

%fitw_m
function [a,b]=fitw_m(xt,yt,w)
% Subfunction for loess_m.m
%
% Original code provided by Rüdiger Röttgers
% Institute for Coastal Research
% Helmholtz Zentrum Geesthacht.
%
% Modified by Jochen Wollschläger
% Institute for Chemistry and Biology of the Marine Environment
% University of Oldenburg
% Contact: jochen.wollschlaeger@uol.de
%
% yt can now also be a matrix (data organized in rows). Allows to fit 
% multiple values at once.
% 
% Changelog:
%
% Version 1.0: Finished 18.01.2019
%--------------------------------------------------------------------------

xt=repmat(xt,size(yt,1),1);
w=repmat(w,size(yt,1),1);

w2 = w.^2;
S = sum(w2,2);
Sx = sum(xt.*w2,2);
Sy = sum(yt.*w2,2);
Sxx = sum(xt.^2.*w2,2);
Sxy = sum(xt.*yt.*w2,2);

delta = S.*Sxx-Sx.^2;
a = (Sxx.*Sy-Sx.*Sxy)./delta;
b = (S.*Sxy-Sx.*Sy)./delta;
end

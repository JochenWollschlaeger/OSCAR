% Script for importing data from the Online Hyperspectral Integrating 
% Cavity Absorption Meter (OSCAR; TriOS GmbH, Germany)
%
% The script imports either the ".dat" or ".csv" files provided by the 
% instrument. It is required that the user uses both the "RAW_DARK" and
% the "RAW_LIGH" files that are obtained after unzipping the .tar file
% downloaded from the OSCAR. The files have to be put in a folder which 
% path can be specified in the "Parameters" section below. 
%
% Furthermore, the script imports information necessary for further
% processing of the data. This includes information about the time of 
% calibrations, reference and sample measurements, as well as the
% coefficients for water absorption, temperature and salinity corrections 
% used in the later processing. Also the diameter of the cavity used can be
% specified. The information have to be provided by the user in the 
% Microsoft Excel file ("User_Input.xlsx"). Instructions for filling out 
% the tables are provided in a seperate sheet in this file.
%
% All data vectors will be interpolated on a common interval (to be 
% determined in the "Parameters" section below) and vector length 
% (determined by the smallest and largest wavelength available). The 
% interpolation is performed by the function 'loess_m.m' that is included at
% the end of the script with its subfunction ('fitw_m.m').The degree of 
% smoothing can be determined separately for the different data in the
% "Parameters" section below.
%
% Finally, the data it will be saved in four textfiles:
%    1. "CalibrationData.txt"
%       Contains the calibration reference and standard measurement as well
%       as the absorption coefficient spectra of the standards used. They 
%       are the mean values of the periods specified in the Excel sheet.
%    2. "Constants.txt"
%       Contains the specified coefficients for water absorption,
%       temperature and salinity correction and cavity radius.
%    3. "ReferenceData.txt"
%       Contains the measurements specified as reference. They are the mean
%       values of the periods specified in the Excel sheet.
%    4. "SampleData.txt"
%       Contains the measurements specified as samples.
% 
% The data in the text files can be edited and corrected, if necessary.
% However, do not change the length of individual vectors, as this causes
% problems for importing the data with the subsequent scripts.
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
%-----------------------------------

clc
clear
BasicPath=pwd;

%% Parameters

%Path of the folder that contains the raw data
RawdataFolder=[BasicPath,'\Rawdata']; %Example, can be replaced by any path

%Name of the folder that will contain the output files
%If not existing, the folder will be created
ResultsFolder=[BasicPath,'\Results']; %Example, can be replaced by any path

%Spectral interval of output data (in nm)
Interval=2;

%Timestamp format in the output files
DesiredFormat='dd.mm.yyyy HH:MM:SS';

%Degree of smoothing for the OSCAR data (loess.m; 0.01-1)
Smoothing_OSCAR=0.05;

%Degree of smoothing for the constants data (loess.m; 0.01-1)
Smoothing_Constants=0.05;

%Degree of smoothing for the absorption coefficient spectra of the standard
%(loess.m; 0.01-1)
Smoothing_Standard=0.075;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FROM HERE, ONLY MAKE MODIFICATIONS IF YOU KNOW WHAT YOU ARE DOING!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check if all required folders and files are available

cd(BasicPath)

%User input file
if ~exist('User_Input.xlsx','file')==2
    error('"User_Input.xlsx" is missing!')
end

%Rawdata
if exist(RawdataFolder,'dir')==7
    cd(RawdataFolder)
    Content=dir;
    Folders=[Content.isdir];
    FileList={Content(~Folders).name}';
    if isempty(FileList)
        error(['No ".csv" or ".dat" files found in the folder ',RawdataFolder,'.'])
    end
    FileType=cellfun(@strsplit,FileList,repmat({'.'},size(FileList)),'UniformOutput',0);
    FileType=vertcat(FileType{:});
    FileType=FileType(:,end);
    FileTypePresent=ismember([{'CSV'} {'DAT'}],FileType);    
    if isequal([true true],FileTypePresent) %both files files present
        error(['Both ".csv" and ".dat" files found. Please have only one type of file in the folder ',RawdataFolder,'.'])      
    end                
else
    error(['The folder "',RawdataFolder,'" is missing!'])
end
clear FileList FileType Content Folders

cd(BasicPath)

%% Import of user input from Excel-file

%Calibration measurements
[~,~,Raw]=xlsread('User_Input.xlsx','Calibrations');
Index=cellfun(@isnan,Raw(:,1),'UniformOutput',0);
Index=cellfun(@any,Index);
Raw=Raw(~Index,:);
Wavelengths_Standard=cell2mat(Raw(1,7:end));
Input=Raw(2:end,:);
if isempty(Input)
    Calib_Name={'No calibration performed'};
    Calib_Ref_Timestamp_Start=datenum('01.01.200100:00','dd.mm.yyyyHH:MM');
    Calib_Ref_Timestamp_End=datenum('01.01.200100:01','dd.mm.yyyyHH:MM');
    Calib_Std_Timestamp_Start=datenum('01.01.200100:02','dd.mm.yyyyHH:MM');
    Calib_Std_Timestamp_End=datenum('01.01.200100:03','dd.mm.yyyyHH:MM');
    AbsCoeff_Standard=NaN(size(Wavelengths_Standard));
else
    Calib_Name=Input(:,1);
    Calib_Ref_Timestamp_Start=datenum(cell2mat([Input(:,2) Input(:,3)]),'dd.mm.yyyyHH:MM');
    Calib_Ref_Timestamp_End=datenum(cell2mat([Input(:,2) Input(:,4)]),'dd.mm.yyyyHH:MM');
    Calib_Std_Timestamp_Start=datenum(cell2mat([Input(:,2) Input(:,5)]),'dd.mm.yyyyHH:MM');
    Calib_Std_Timestamp_End=datenum(cell2mat([Input(:,2) Input(:,6)]),'dd.mm.yyyyHH:MM');
    AbsCoeff_Standard=cell2mat(Input(:,7:end));
end
clear Raw Index Input

%Reference measurements
[~,~,Raw]=xlsread('User_Input.xlsx','Reference measurements');
Index=cellfun(@isnan,Raw(:,1),'UniformOutput',0);
Index=cellfun(@any,Index);
Raw=Raw(~Index,:);
Raw=Raw(2:end,:);
if isempty(Raw)
    error('No reference measurements specified!')
else
    Ref_Name=Raw(:,1);
    Ref_Timestamp_Start=datenum(cell2mat([Raw(:,2) Raw(:,3)]),'dd.mm.yyyyHH:MM');
    Ref_Timestamp_End=datenum(cell2mat([Raw(:,2) Raw(:,4)]),'dd.mm.yyyyHH:MM');
end
clear Raw Index

%Sample measurements
[~,~,Raw]=xlsread('User_Input.xlsx','Sample measurements');
Index=cellfun(@isnan,Raw(:,1),'UniformOutput',0);
Index=cellfun(@any,Index);
Raw=Raw(~Index,:);
Raw=Raw(2:end,:);
if isempty(Raw)
    error('No sample measurements specified!')
else
    Sample_Name=Raw(:,1);
    Sample_Timestamp_Start=datenum(cell2mat([Raw(:,2) Raw(:,3)]),'dd.mm.yyyyHH:MM');
    Sample_Timestamp_End=datenum(cell2mat([Raw(:,4) Raw(:,5)]),'dd.mm.yyyyHH:MM');
end
clear Raw Index

%Constants
[~,~,Raw]=xlsread('User_Input.xlsx','Constants');
Index=cellfun(@isnan,Raw(:,1),'UniformOutput',0);
Index=cellfun(@any,Index);
Raw=Raw(~Index,:);
Raw=cell2mat(Raw(2:end,:));
if isempty(Raw)
    error('No constants data given!')
end
if isnan(Raw(:,1))
    error('No wavelengths available for the constants!')
else
    ConstantsWavelengths=transpose(Raw(:,1));
end
if isnan(Raw(:,2))
    error('No water absorption coeffients given!')
else
    AbsCoeff_Water=transpose(Raw(:,2));
end
if isnan(Raw(:,3))
    error('No temperature correction coefficients given!')
else
    TempCoeff=transpose(Raw(:,3));
end
if isnan(Raw(:,4))
    error('No salinity correction coefficients given!')
else
    SaltCoeff=transpose(Raw(:,4));
end
if isnan(Raw(:,5))
    error('The temperature at which the absorption coefficients of pure water has been determined is missing!')
else
    AbsCoeff_WaterT=transpose(Raw(1,5));
end
if isnan(Raw(1,6))
    error('The radius of the cavity is not specified!')
else
    r=Raw(1,6);
end
if isnan(Raw(1,7))
    error('The radius of the cavity minus the radius of the quartz ball is not specified!')
else
    r0=Raw(1,7);
end
clear Raw Index

%% Checking the timestamps given in the manual input

%Start-time of measurement smaller than end-time of measurement?
if sum(~(Calib_Ref_Timestamp_Start<Calib_Ref_Timestamp_End))>0
    error('Check Calibration Reference times in Excel file');
end
if sum(~(Calib_Std_Timestamp_Start<Calib_Std_Timestamp_End))>0
    error('Check Calibration Standard times in Excel file');
end
if sum(~(Ref_Timestamp_Start<Ref_Timestamp_End))>0
    error('Check Reference times in Excel file');
end
if sum(~(Sample_Timestamp_Start<Sample_Timestamp_End))>0
    error('Check Calibration Reference times in Excel file');
end

%% Import OSCAR data

cd(RawdataFolder)

%Creating a list of all files present in the folder containing the raw data
Content=dir;
Folders=[Content.isdir];
FileList={Content(~Folders).name}';
clear Content Folders

%Opening the files one after another and merge the data
OSCAR_Timestamp=cell(length(FileList),1);
OSCAR_SpectrumType=cell(length(FileList),1);
OSCAR_IntTime=cell(length(FileList),1);
OSCAR_Wavelengths=cell(length(FileList),1);
OSCAR_Data=cell(length(FileList),1);

if isequal([true false],FileTypePresent) % Data is given as .CSV files
    for i=1:size(FileList,1)
        
        fid=fopen(FileList{i});
        Format=repmat({'%s'},1,500);
        Header=textscan(fid,cell2mat(Format),1,'delimiter',',','EmptyValue',NaN);
        Header=cellfun(@strrep,Header,repmat({''''},size(Header)),repmat({''},size(Header)));
        Header=str2double(Header);
        Index=isnan(Header);
        Format(~Index)={'%f'};
        FileContent=textscan(fid,cell2mat(Format),'delimiter',',','EmptyValue',NaN);
        
        %Obtain the timestamps
        OSCAR_Timestamp{i,1}=FileContent{1};
        
        %Obtain wavelengths and measured spectra
        OSCAR_Wavelengths{i,1}=Header(~Index);
        OSCAR_Data{i,1}=cell2mat(FileContent(~Index));
        fclose(fid);
    end
    clear i Index fid ans
    OSCAR_Wavelengths=unique(cell2mat(OSCAR_Wavelengths),'rows');
    OSCAR_Timestamp=vertcat(OSCAR_Timestamp{:});
    OSCAR_Timestamp=datenum(OSCAR_Timestamp,'yyyy/mm/dd HH:MM:SS');
    OSCAR_Data=cell2mat(OSCAR_Data);
    
elseif isequal([false true],FileTypePresent) %Data is given as .DAT files    
    for i=1:size(FileList,1)
        fid=fopen(FileList{i});
        RawData=textscan(fid,'%s%s%s%s%*[^\n]','delimiter',' ','CommentStyle','[','EmptyValue',NaN);
        fclose(fid);
        clear fid ans

        %Obtain the timestamps
        Index=ismember(RawData{1}, {'DateTime'});
        OSCAR_Timestamp{i,1}=datenum(cell2mat([RawData{3}(Index) RawData{4}(Index)]),'yyyy-mm-ddHH:MM:SS');
        clear Index
        
        %Obtain the spectrum type
        Index=ismember(RawData{1}, {'SpectrumType'});
        OSCAR_SpectrumType{i,1}=RawData{3}(Index);        
        clear Index
        
        %Obtain the integration time
        Index=ismember(RawData{1}, {'IntegrationTime'});
        OSCAR_IntTime{i,1}=str2double(RawData{3}(Index));        
        clear Index

        %Obtain wavelengths and measured spectra
        Wavelengths=str2double(RawData{1});
        Data=str2double(RawData{2});
        Index=isnan(Wavelengths);
        Wavelengths=Wavelengths(~Index);
        Data=Data(~Index);
        Index=Wavelengths==0;
        Wavelengths=Wavelengths(~Index);
        Data=Data(~Index);
        WavelengthsUnique=unique(Wavelengths,'stable');
        Spectra=NaN(length(OSCAR_Timestamp{i,1}),length(WavelengthsUnique));
        for j=1:length(WavelengthsUnique)
            Index=Wavelengths==WavelengthsUnique(j);
            OneWavelength=Data(Index);
            Spectra(:,j)=OneWavelength;
        end
        OSCAR_Wavelengths{i,1}=WavelengthsUnique';
        OSCAR_Data{i,1}=Spectra;
        clear j WavelengthsUnique Index RawData Data OneWavelength Spectra...
              Wavelengths     
    end
    OSCAR_Wavelengths=unique(cell2mat(OSCAR_Wavelengths),'rows');
    OSCAR_SpectrumType=vertcat(OSCAR_SpectrumType{:});
    OSCAR_IntTime=cell2mat(OSCAR_IntTime);
    OSCAR_Timestamp=cell2mat(OSCAR_Timestamp);
    OSCAR_Data=cell2mat(OSCAR_Data);  
end
clear i RawData fid ans FileList FileTypePresent

%% Calculating the measurement spectrum

%Normalizing all data on integration time
OSCAR_Data=OSCAR_Data./repmat(OSCAR_IntTime,1,size(OSCAR_Data,2));
clear OSCAR_IntTime

%Division in dark and light spectra
Index=ismember(OSCAR_SpectrumType,{'RAW_Dark'});
OSCAR_Timestamp_Dark=OSCAR_Timestamp(Index,:);
OSCAR_Timestamp_Light=OSCAR_Timestamp(Index,:);
OSCAR_Data_Dark=OSCAR_Data(Index,:);
OSCAR_Data_Light=OSCAR_Data(~Index,:);
clear Index OSCAR_SpectrumType

%Check whether each light spectrum has a corresponding dark spectrum
if size(OSCAR_Data_Light,1)==size(OSCAR_Data_Dark,1)
    if sum(OSCAR_Timestamp_Dark-OSCAR_Timestamp_Light)==0
    else
        error('There are differences in the timestamps of dark and light spectra!')
    end
else
    error('Something is wrong with the number of dark and light spectra!')
end

%Subtracting dark spectra from light spectra
OSCAR_Timestamp=OSCAR_Timestamp_Light;
OSCAR_Data=OSCAR_Data_Light-OSCAR_Data_Dark;
clear OSCAR_Data_Dark OSCAR_Data_Light OSCAR_Timestamp_Dark OSCAR_Timestamp_Light

%Deletion of doubled timestamps and corresponding data
[OSCAR_Timestamp,Index,~]=unique(OSCAR_Timestamp);
OSCAR_Data=OSCAR_Data(Index,:);
clear Index

%% Restricting the data to wavelengths that are available for all parameters

%Determination of minimum/maximum wavelength available
WavelengthMinimum=max([min(ConstantsWavelengths),min(Wavelengths_Standard),min(OSCAR_Wavelengths)]);
WavelengthMaximum=min([max(ConstantsWavelengths),max(Wavelengths_Standard),max(OSCAR_Wavelengths)]);

%Restricting the "Constants"-data
Index=ConstantsWavelengths>=WavelengthMinimum&ConstantsWavelengths<=WavelengthMaximum;
ConstantsWavelengths=ConstantsWavelengths(:,Index);
AbsCoeff_Water=AbsCoeff_Water(:,Index);
TempCoeff=TempCoeff(:,Index);
SaltCoeff=SaltCoeff(:,Index);
clear Index

%Restricting the absorption coefficient spectra of the calibration standards
Index=Wavelengths_Standard>=WavelengthMinimum&Wavelengths_Standard<=WavelengthMaximum;
Wavelengths_Standard=Wavelengths_Standard(:,Index);
AbsCoeff_Standard=AbsCoeff_Standard(:,Index);
clear Index

%Restricting the OSCAR data
Index=OSCAR_Wavelengths>=WavelengthMinimum&OSCAR_Wavelengths<=WavelengthMaximum;
OSCAR_Wavelengths=OSCAR_Wavelengths(:,Index);
OSCAR_Data=OSCAR_Data(:,Index);
clear Index

%% Interpolation of the spectra according to the parameters set

%Setting wavelength range and interval
Wavelengths=WavelengthMinimum:Interval:WavelengthMaximum;
clear WavelengthMinimum Interval WavelengthMaximum

%Interpolating the "Constants"-data
AbsCoeff_Water=loess_m(ConstantsWavelengths,AbsCoeff_Water,Wavelengths,Smoothing_Constants);
TempCoeff=loess_m(ConstantsWavelengths,TempCoeff,Wavelengths,Smoothing_Constants);
SaltCoeff=loess_m(ConstantsWavelengths,SaltCoeff,Wavelengths,Smoothing_Constants);
clear ConstantsWavelengths alpha_Constants

%Interpolating the absorption coefficient spectra of the calibration standards
New_AbsCoeff_Standard=NaN(size(AbsCoeff_Standard,1),size(Wavelengths,2));
for i=1:size(AbsCoeff_Standard,1)
    New_AbsCoeff_Standard(i,:)=loess_m(Wavelengths_Standard,AbsCoeff_Standard(i,:),Wavelengths,Smoothing_Standard);
end
AbsCoeff_Standard=New_AbsCoeff_Standard;
clear i New_AbsCoeff_Standard Wavelengths_Standard alpha_Standard

%Interpolating the OSCAR data
New_OSCAR_Data=zeros(size(OSCAR_Data,1),size(Wavelengths,2));
InterpolationsToDo=size(OSCAR_Data,1);
for i=1:size(OSCAR_Data,1)
    New_OSCAR_Data(i,:)=loess_m(OSCAR_Wavelengths,OSCAR_Data(i,:),Wavelengths,Smoothing_OSCAR);
    InterpolationsToDo=InterpolationsToDo-1;
    display(InterpolationsToDo)
end
OSCAR_Data=New_OSCAR_Data;
clear i New_OSCAR_Data InterpolationsToDo OSCAR_Wavelengths alpha_OSCAR

cd(BasicPath)

%% Division of the data into the different types of measurement

%Calibration Reference measurements
Calib_Ref_Data=NaN(size(Calib_Name,1),size(OSCAR_Data,2));
Calib_Ref_N=NaN(size(Calib_Name,1),1);
for i=1:size(Calib_Name,1)
    Index=OSCAR_Timestamp>Calib_Ref_Timestamp_Start(i)&OSCAR_Timestamp<Calib_Ref_Timestamp_End(i);
    Calib_Ref_Data(i,:)=mean(OSCAR_Data(Index,:),1,'omitnan');
    Calib_Ref_N(i,:)=numel(OSCAR_Data(Index,1));
end
clear i Index

%Calibration Standard measurements
Calib_Standard_Data=NaN(size(Calib_Name,1),size(OSCAR_Data,2));
Calib_Standard_N=NaN(size(Calib_Name,1),1);
for i=1:size(Calib_Name,1)
    Index=OSCAR_Timestamp>Calib_Std_Timestamp_Start(i)&OSCAR_Timestamp<Calib_Std_Timestamp_End(i);
    Calib_Standard_Data(i,:)=mean(OSCAR_Data(Index,:),1,'omitnan');
    Calib_Standard_N(i,:)=numel(OSCAR_Data(Index,1));
end
clear i Index

%Reference measurements
Ref_Data=NaN(size(Ref_Name,1),size(OSCAR_Data,2));
Ref_N=NaN(size(Ref_Name,1),1);
for i=1:size(Ref_Name,1)
    Index=OSCAR_Timestamp>Ref_Timestamp_Start(i)&OSCAR_Timestamp<Ref_Timestamp_End(i);
    Ref_Data(i,:)=mean(OSCAR_Data(Index,:),1,'omitnan');
    Ref_N(i,:)=numel(OSCAR_Data(Index,1));
end
clear i Index

%Sample measurements
Sample_Data=cell(size(Sample_Name,1),1);
Sample_Timestamp=cell(size(Sample_Name,1),1);
Sample_Name_extended=cell(size(Sample_Name,1),1);
for i=1:size(Sample_Name,1)
    Index=OSCAR_Timestamp>Sample_Timestamp_Start(i)&OSCAR_Timestamp<Sample_Timestamp_End(i);
    Sample_Data{i,1}=OSCAR_Data(Index,:);
    Sample_Timestamp{i,1}=OSCAR_Timestamp(Index,:);
    Sample_Name_extended{i,1}=repmat(Sample_Name(i,1),sum(Index),1);
end
clear i Index OSCAR_Data OSCAR_Timestamp

Sample_Data=cell2mat(Sample_Data);
Sample_Timestamp=cell2mat(Sample_Timestamp);
Sample_Name=vertcat(Sample_Name_extended{:});
clear Sample_Name_extended

%% Converting the timestamps back in string format

Calib_Ref_Timestamp_Start=cellstr(datestr(Calib_Ref_Timestamp_Start,DesiredFormat));
Calib_Ref_Timestamp_End=cellstr(datestr(Calib_Ref_Timestamp_End,DesiredFormat));

Calib_Std_Timestamp_Start=cellstr(datestr(Calib_Std_Timestamp_Start,DesiredFormat));
Calib_Std_Timestamp_End=cellstr(datestr(Calib_Std_Timestamp_End,DesiredFormat));

Ref_Timestamp_Start=cellstr(datestr(Ref_Timestamp_Start,DesiredFormat));
Ref_Timestamp_End=cellstr(datestr(Ref_Timestamp_End,DesiredFormat));

Sample_Timestamp=cellstr(datestr(Sample_Timestamp,DesiredFormat));

%% Saving the data

%Going to or creating the folder where the working data will be saved
if exist(ResultsFolder,'dir')==0
    mkdir(ResultsFolder)
end

cd(ResultsFolder)

%Constants data
SaveFileName='Constants.txt';

HeaderFormat=['%s',repmat('\t%d',1,size(Wavelengths,2)),'\r\n'];
Header=[{'Wavelengths [nm]'} num2cell(Wavelengths)];

DataFormat=['%s',repmat('\t%f',1,size(Wavelengths,2)),'\r\n'];
Data=[[{'Absorption coefficient water [m-1]'} num2cell(AbsCoeff_Water)];...
      [{'Temperature coefficient [m-1]'} num2cell(TempCoeff)];...
      [{'Salinity coefficient [m-1]'} num2cell(SaltCoeff)]];

fid=fopen(SaveFileName,'w');
fprintf(fid,'%s\r\n',['Radius cavity [m]: ',num2str(r)]);
fprintf(fid,'%s\r\n',['Radius cavity - radius light source [m]: ',num2str(r0)]);
fprintf(fid,'%s\r\n',['Temperature of water absorption coefficient spectrum [°C]: ',num2str(AbsCoeff_WaterT)]);
fprintf(fid,'\r\n');
fprintf(fid,HeaderFormat,Header{1,:});
for i=1:size(Data,1)
    fprintf(fid,DataFormat,Data{i,:});
end
fclose(fid);
clear i ans fid SaveFileName Header HeaderFormat Data DataFormat

%Calibration data
SaveFileName='CalibrationData.txt';

HeaderFormat=['%s\t%s\t%s\t%s\t%s\t%s',repmat('\t%d',1,size(Wavelengths,2)),'\r\n'];
Header=[{'Calibration No./Name'} {'Label [DO NOT RENAME!]'} {['Timestamp Start [',DesiredFormat,']']}...
        {['Timestamp End [',DesiredFormat,']']} {'Temperature [°C]'} {'N'} num2cell(Wavelengths)];

DataFormat=['%s\t%s\t%s\t%s\t%.2f\t%.0f',repmat('\t%f',1,size(Wavelengths,2)),'\r\n'];
Calib_Ref_Data=[Calib_Name repmat({'Calib_Ref'},size(Calib_Name)) Calib_Ref_Timestamp_Start Calib_Ref_Timestamp_End...
                num2cell(NaN(size(Calib_Name))) num2cell(Calib_Ref_N) num2cell(Calib_Ref_Data)];
Calib_Standard_Data=[Calib_Name repmat({'Calib_Standard'},size(Calib_Name)) Calib_Std_Timestamp_Start...
                     Calib_Std_Timestamp_End num2cell(NaN(size(Calib_Name))) num2cell(Calib_Standard_N) num2cell(Calib_Standard_Data)];
AbsCoeff_Standard_Data=[Calib_Name repmat({'AbsCoeff_Standard'},size(Calib_Name)) repmat({'N.A.'},size(Calib_Name)) repmat({'N.A.'},size(Calib_Name)) num2cell(NaN(size(Calib_Name))) num2cell(NaN(size(Calib_Name))) num2cell(AbsCoeff_Standard)];
Data=[Calib_Ref_Data;Calib_Standard_Data;AbsCoeff_Standard_Data];

fid=fopen(SaveFileName,'w');
fprintf(fid,HeaderFormat,Header{1,:});
for i=1:size(Data,1)
    fprintf(fid,DataFormat,Data{i,:});
end
fclose(fid);
clear i fid ans SaveFileName Header HeaderFormat Data DataFormat...
      Calib_Ref_Data Calib_Standard_Data AbsCoeff_Standard_Data

%Reference data
SaveFileName='ReferenceData.txt';

HeaderFormat=['%s\t%s\t%s\t%s\t%s',repmat('\t%d',1,size(Wavelengths,2)),'\r\n'];
Header=[{'Reference No./Name'} {['Timestamp Start [',DesiredFormat,']']} {['Timestamp End [',DesiredFormat,']']}...
        {'Temperature [°C]'} {'N'} num2cell(Wavelengths)];

DataFormat=['%s\t%s\t%s\t%.2f\t%.0f',repmat('\t%f',1,size(Wavelengths,2)),'\r\n'];
Data=[Ref_Name Ref_Timestamp_Start Ref_Timestamp_End num2cell(NaN(size(Ref_Name))) num2cell(Ref_N) num2cell(Ref_Data)];

fid=fopen(SaveFileName,'w');
fprintf(fid,HeaderFormat,Header{1,:});
for i=1:size(Data,1)
    fprintf(fid,DataFormat,Data{i,:});
end
fclose(fid);
clear i fid ans SaveFileName Header HeaderFormat Data DataFormat
  
%Sample data
SaveFileName='SampleData.txt';

HeaderFormat=['%s\t%s\t%s\t%s',repmat('\t%d',1,size(Wavelengths,2)),'\r\n'];
Header=[{'Sample No./Name'} {['Timestamp [',DesiredFormat,']']} {'Temperature [°C]'} {'Salinity'} num2cell(Wavelengths)];

DataFormat=['%s\t%s\t%.2f\t%.0f',repmat('\t%f',1,size(Wavelengths,2)),'\r\n'];
Data=[Sample_Name Sample_Timestamp num2cell(NaN(size(Sample_Name))) num2cell(NaN(size(Sample_Name))) num2cell(Sample_Data)];

fid=fopen(SaveFileName,'w');
fprintf(fid,HeaderFormat,Header{1,:});
for i=1:size(Data,1)
    fprintf(fid,DataFormat,Data{i,:});
end
fclose(fid);
clear i ans fid SaveFileName Header HeaderFormat Data DataFormat

cd(BasicPath)

Properties=struct('WindowStyle','non-modal','Interpreter','tex');
msgbox('\bf \fontsize{12} \fontname{Arial} Script finished',Properties)
clear Properties

%% Non build-in functions required from the script

%loess
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

%fitw
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
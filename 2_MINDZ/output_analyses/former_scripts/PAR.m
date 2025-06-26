% This script compute PAR profiles according to surface irradiance and
% dissolved coloured matter concentrations, approximated by a relationship
% with surface salinity in the GSL.
%
% MapsF 2012

%% Import data


% Surface irradiance

inPAR = input('Input PAR file name? ');

ncid = netcdf.open(inPAR,'nowrite');


xdimid = netcdf.inqDimID(ncid,'lon');

ydimid = netcdf.inqDimID(ncid,'lat');

tdimid = netcdf.inqDimID(ncid,'time');

[~, xlen] = netcdf.inqDim(ncid,xdimid);

[~, ylen] = netcdf.inqDim(ncid,ydimid);

[~, tlen] = netcdf.inqDim(ncid,tdimid);


lonid = netcdf.inqVarID(ncid,'lon');

latid = netcdf.inqVarID(ncid,'lat');

tid   = netcdf.inqVarID(ncid,'time');

% Variable vbdsf (Visible Beam Downstream Solar Flux; W.m^-2) is compressed
Iid   = netcdf.inqVarID(ncid,'vbdsf');

offset = double( netcdf.getAtt(ncid,Iid,'add_offset') );

scale  = double( netcdf.getAtt(ncid,Iid,'scale_factor') );


lon  = double( netcdf.getVar(ncid,lonid) );

lat  = double( netcdf.getVar(ncid,latid) );

time = double( netcdf.getVar(ncid,tid) );


start = [0 0 0];

count = [xlen ylen tlen];

I = double( netcdf.getVar(ncid,Iid,start,count) )*scale+offset;


% Weighted mean depth & surface (10m) salinity

inkrill = input('Input file name WMD vs Salinity : ');

newdata = importdata(inkrill);

% Create new variables from those fields
vars = fieldnames(newdata); % 'data' & 'textdata'
for i = 1:length(vars)
    assignin('base', vars{i}, newdata.(vars{i}));
end

cil = data(:,1); % Depth of the bottom of the CIL

S   = data(:,2); % Surface salinity

Tmn = data(:,4); % Temp. at M. norvegica WMD

Ttr = data(:,5); % Temp. at T. raschii WMD

WMD = data(:,6:7);  % WMD M. norvegica & T. raschii

clear newdata vars data textdata


%% Parameters

kw = 0.04; % Constant attenuation coefficient of sea water (m^-1)

W2microEin = 4.5e-6; % Conversion factor between W.m^-2 and µEin.s^-1.m^-2


%% Compute PAR profile

Iave = mean(mean(mean(I(1:4,1:2,491:4:615),3))); % Average July light intensity at surface

kp = -0.0364.*S+1.1942; % Variable attenuation coefficient of colored dissolved matter (m^-1)

PARmix = 0.45*Iave.*exp(-(kw+kp).*10); % PAR under the surface mixed layer (10m) where kp counts

PARwmd = [PARmix PARmix].*exp(-kw.*(WMD-10)); % PAR at the WMD, from the constant attenuation of sea water below kp influance

PARein = PARwmd*W2microEin; % Convert in µEin.s^-1.m^-2

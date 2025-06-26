% krill_load
%
% Load data for the spatial comparison between FTLE and observed krill density in the GSL.
%
% MapsF 2013

%% Load data

% Krill file

if ~exist('krillfile','var')

    krillfile = ['~/Documents/JOB/NEMO/Data/' input('\nInput krill file name (.csv): ')];

    datain = importdata(krillfile);

    % Create new variables in the base workspace from those fields
    for i = 1:size(datain.colheaders,2)
        assignin('base', genvarname(datain.colheaders{i}), datain.data(:,i));
    end

    % Chose main variable
    disp(' ')
    disp(datain.colheaders)

    i = input('Which of the variable above is the variable of interest (give index #)? ');

    year = datain.data(:,1);

    data = datain.data(:,i);

    clear datain

end


% Model coordinates file: lon/lat from krill files converted in model x/y
% with routine lonlat_ij.f90

if ~exist('xy','var')

    xy = load('~/Documents/JOB/NEMO/Data/xy.txt');

    x = xy(:,1); y = xy(:,2);

end


% FTLE data
disp(' ')
trfile = input('Input T. raschii FTLE file name (without ''.mat''): ');

load(trfile)

Tr = FTLE;

disp(' ')
sfile = input('Input surface tracers FTLE file name (without ''.mat''): ');

load(sfile)

Surf = FTLE; clear FTLE;


% Load the matrices of distances between the dense aggregations and the
% closest FTLE patches (Computed in R with the script FTLE_Matlab_dist.R)

load('~/Documents/JOB/OUTPUT/POST/NEMO_Gulf_Slope_distance.mat')

gulf_dist = result;

load('~/Documents/JOB/OUTPUT/POST/NEMO_Anticosti_Slope_distance.mat')

anti_dist = result;

load([trfile '_distance'])

Tr_dist = result;

load([sfile '_distance'])

Surf_dist = result;

clear result;


% Surface currents (every days of August 2006-2010)
ncid = netcdf.open('~/Documents/JOB/FORCING/l_2011_m2_ave_UV.nc','nowrite');

uid = netcdf.inqVarID(ncid,'U');

vid = netcdf.inqVarID(ncid,'V');

start = [0 0 0 0]; count = [197 234 1 155]; % Extract surface layer

U = squeeze(double( netcdf.getVar(ncid,uid,start,count) )); U(U==0) = nan;

V = squeeze(double( netcdf.getVar(ncid,vid,start,count) )); V(V==0) = nan;


% Lat/Lon of partikrill trajectories

if ~exist('intraj','var')

    disp(' ')
    intraj = input('Enter input filename for backward particles trajectories: ');

    if ~isempty(intraj)

        ncid = netcdf.open(['~/Documents/JOB/OUTPUT/' intraj],'nowrite');

        lonid = netcdf.inqVarID(ncid,'lon');

        latid = netcdf.inqVarID(ncid,'lat');

        tid = netcdf.inqVarID(ncid,'time');

        start = [0 0 31 31]; count = [197 234 15 30];

        llon = double( netcdf.getVar(ncid,lonid,start,count) ); llon(llon>1e9) = nan;

        llat = double( netcdf.getVar(ncid,latid,start,count) ); llat(llat>1e9) = nan;

        lltime = double( netcdf.getVar(ncid,tid,31,15) );
    end

end


% Model bathymetry

if ~exist('slope','var')

    load('~/Documents/JOB/NEMO/Data/NEMO_bottom_depth.mat')

    slope(isnan(Tr(:,:,1))) = nan;
    bottom(bottom==0) = nan;

    bott = nan(size(bottom));

    
% Coastline
    load ~/Documents/JOB/NEMO/Data/NEMO_coastline.mat

end

return

% krill_spatial_comp
%
% Compute spatial comparison between FTLE and observed krill density in the GSL.
%
% MapsF 2013


%% Prepare data for spatial analysis
%  Data has been loaded beforehand
%
%  year  = year of observations
%  data  = observed accoustic data
%  x     = coordinate of the data point on the model grid
%  y     = ""
%  Lon   = Longitude at the center of each cell from the model grid
%  Lat   = Latitude at the center of each cell from the model grid
%  Tr    = FTLE results for particles with a T. raschii DVM
%  Surf  = FTLE results for particules distributed around 15m depth
%  Slope = slope of the model grid

% Select period
mi = true(size(data));

% Select region
ri = x>-999 & y>-999 & Lon>-70 & Lat>46;

% Get rid of NaNs in krill data
ki = isfinite(data);

% Create explicative variables (FTLE)
ftle1 = Tr; ftle2 = Surf;

for i = 1:size(ftle1,3)
    ftle1(:,:,i) = ftle1(:,:,i)/max(max(ftle1(:,:,i)));
end

for i = 1:size(ftle2,3)
    ftle2(:,:,i) = ftle2(:,:,i)/max(max(ftle2(:,:,i)));
end

ftle1 = nanmean(ftle1,3);

ftle2 = nanmean(ftle2,3);

% Get FTLE and bathymetry values at the locations of observations
dumm1 = nan(size(x)); dumm2 = dumm1; dummz = dumm1;

for i = 1:length(x)

    dumm1(i,1) = ftle1(max(1,ceil(x(i))),max(1,ceil(y(i))));
    dumm2(i,1) = ftle2(max(1,ceil(x(i))),max(1,ceil(y(i))));

    dummz(i,1) = slope(max(1,ceil(x(i))),max(1,ceil(y(i))));

end

fi = isfinite(dumm1);

% Select data

Lat1 = Lat(mi&ri&ki&fi);

Lon1 = Lon(mi&ri&ki&fi);

Year1 = year(mi&ri&ki&fi);

dummy = data(mi&ri&ki&fi);
% Assumes a dB to biomass conversion factor of -70 dB g^-1 (McQuinn et al. 2013)
% ( 1 dB = 10*log10(x) )
dummy = dummy./10^(-70/10); % g m^-2


n = length(dummy);

%dummx(~mi|~ri|~ki|~fi) = [];
dumm1(~mi|~ri|~ki|~fi) = [];
dumm2(~mi|~ri|~ki|~fi) = [];

dummz(~mi|~ri|~ki|~fi) = [];


%% Compute correspondence for predicted (FTLE) and observed (krill abundance) aggregation areas

q1 = zeros(10,1); q2 = q1; q3 = q1; qz = q1;

dx = qz; ds = dx; dz = dx; 

PX11 = dz; PX = PX11; PS11 = PX; PS = PX; PZ11 = PX; PZ = PX;

bin1 = zeros(n,10); bin2 = bin1; binz = bin1;

mask = isfinite(ftle1)&lon>-70;

% Extract upper quantile of zooplankton density = dense aggregations
%qy = quantile(dummy(dummy>0),0.95);

% Threshold for dense aggregations is 100 g m^-3 accroding to the
% bioenergetics model of blue whales from Goldbogen et al. (2011)
% Then our threshold is 100 g m^-2, i.e. the minimum vertically integrated
% abundance required to produce such a density by any patch forming mechanism...
qy = 100; % g m^-2
biny = dummy; biny(biny<qy) = 0; biny(biny>=qy) = 1;

for i = 1:10

    % Extract upper quantile of FTLE values = convergence areas
    q1(i) = quantile(reshape(ftle1(mask),numel(ftle1(mask)),1),0.45+i*0.05);
    q2(i) = quantile(reshape(ftle2(mask),numel(ftle2(mask)),1),0.45+i*0.05);

    % Extract upper quantile of slope = bathymetric discontinuity
    qz(i) = quantile(reshape(slope(mask),numel(slope(mask)),1),0.45+i*0.05);

    % Compute binary values
    bin1(:,i) = dumm1; bin1(dumm1<q1(i),i) = 0; bin1(dumm1>=q1(i),i) = 1;

    bin2(:,i) = dumm2; bin2(dumm2<q2(i),i) = 0; bin2(dumm2>=q2(i),i) = 1;

    binz(:,i) = dummz; binz(dummz<qz(i),i) = 0; binz(dummz>=qz(i),i) = 1;


    % Compute ratio of observed match probability to random match probability
    P1 = sum(biny==1)/n; % Probability to observe an aggregation

    PX11(i) = sum(bin1(:,i)==1&biny==1)/sum(biny==1); % P that an observed aggregation matches a zone of DVM convergence

    PX(i) =  PX11(i)/(P1*(1-0.45-i*0.05)); % 1-0.45-i*0.05 = P that a cell from the model is a zone of DVM convergence

    PS11(i) = sum(bin2(:,i)==1&biny==1)/sum(biny==1); % P that an observed aggregation matches a zone of surface convergence

    PS(i) =  PS11(i)/(P1*(1-0.45-i*0.05)); % 1-0.45-i*0.05 = P that a cell from the model is a zone of surface convergence

    PZ11(i) = sum(binz(:,i)==1&biny==1)/sum(biny==1); % P that an observed aggregation matches a zone of bathymetric discontinuity

    PZ(i) =  PZ11(i)/(P1*(1-0.45-i*0.05));


    % Compute mean density observed inside convergence zones
    dx(i) = mean(dummy(dummy>0&bin1(:,i)==1)); % with DVM

    ds(i) = mean(dummy(dummy>0&bin2(:,i)==1)); % surface  

    % Compute mean density observed inside zones of bathymetric
    % discontinuity
    dz(i) = mean(dummy(dummy>0&binz(:,i)==1));

end

return

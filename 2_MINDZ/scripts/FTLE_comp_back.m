% This routine computes Finite-Time Lyapunov Exponents (FTLE) fields from a
% matrix of lagrangian particles positions evolving BACKWARD in time.
%
% http://mmae.iit.edu/shadden/LCS-tutorial/computation.html
%
% MapsF 2012

function ret = FTLE_comp_back(infile,dt)

%tic

%% Load input NetCDF files

ret = -99;

%infile = input('Input file name ? ');

%dt = input('Time step for computing FTLE (in days, even)? ');


ncid = netcdf.open(['/home/MapsF/Documents/JOB/OUTPUT/' infile],'nowrite');


pdimid = netcdf.inqDimID(ncid,'part');

tdimid = netcdf.inqDimID(ncid,'time');

[~, plen] = netcdf.inqDim(ncid,pdimid);

[~, tlen] = netcdf.inqDim(ncid,tdimid);


lonid = netcdf.inqVarID(ncid,'lon');

latid = netcdf.inqVarID(ncid,'lat');

dxid  = netcdf.inqVarID(ncid,'dx');

dyid  = netcdf.inqVarID(ncid,'dy');

xid   = netcdf.inqVarID(ncid,'xt');

yid   = netcdf.inqVarID(ncid,'yt');

tid   = netcdf.inqVarID(ncid,'time');


time = double( netcdf.getVar(ncid,tid) );


m  = 197; n  = 234; % # of cells in the model grid
%mi = 0; ni = 0;
%mf = m; nf = n;

% Gulf of St Lawrence
mi = 24;  ni = 120;
mf = 190; nf = 227;

% Around Anticosti
%mi = 74;  ni = 169; % Initial coordinates for the (sub)sampling grid for computing FTLEs
%mf = 145; nf = 213; % Final coordinates for the (sub)sampling grid

% Gulf of Maine
%mi = 5;  ni = 40;
%mf = 95; nf = 125;


start = [0 0];

count = [m n];

lon  = double( netcdf.getVar(ncid,lonid,start,count) );

lat  = double( netcdf.getVar(ncid,latid,start,count) );

dx   = double( netcdf.getVar(ncid,dxid,start,count) );

dy   = double( netcdf.getVar(ncid,dyid,start,count) );


start = [mi ni 0 0];

count = [mf-mi nf-ni plen tlen];

x = double( netcdf.getVar(ncid,xid,start,count) ); x(x>1e9) = nan;

y = double( netcdf.getVar(ncid,yid,start,count) ); y(y>1e9) = nan;


%% Get the correct timesteps & convert coordinates in distances
% Data format: X|Y(i,j,T0,T)

X = nan(size(x,1),size(x,2),size(x,3),2); Y = X;

for i = 1:size(x,3)

    X(:,:,i,1) = x(:,:,i,i);

    X(:,:,i,2) = x(:,:,i,i+round(dt*0.5));

    Y(:,:,i,1) = y(:,:,i,i);

    Y(:,:,i,2) = y(:,:,i,i+round(dt*0.5));

end

clear x y

X(X<1) = 1; Y(Y<1) = 1;

X(X>m) = m; Y(Y>n) = n;

% Interpolate position in km according to the dx & dy matrices

DX = zeros(size(dx,1)+1, size(dx,2));

DX(2:end,:) = cumsum(dx,1);


DY = zeros(size(dy,1), size(dy,2)+1);

DY(:,2:end) = cumsum(dy,2);

% Select the area inside the Gulf of St Lawrence based on model's grid incdices
a = (32-7)/(146-118); b = 7-a*118;

% Select the area around Anticosti Island
% %a = (160-210)/(150-60); b = 210-a*60;
%a = (-9-41)/(76+14); b = 41-a*(-14);

for i = 1:size(X,1)
    j = floor(a*i+b);
    X(i,1:j,1,1) = nan;
end

[i,j] = find( isfinite( squeeze(X(:,:,1,1)) ) );

XX = nan(size(X));

YY = nan(size(Y));

for k = 1:length(i)

    for l = 1:size(X,3)*size(X,4)

            XX(i(k),j(k),l) = DX( ceil(X(i(k),j(k),l)), ceil(Y(i(k),j(k),l)) )...
                             + dx( ceil(X(i(k),j(k),l)), ceil(Y(i(k),j(k),l)) )...
                             .*(X(i(k),j(k),l)-floor(X(i(k),j(k),l)));

            YY(i(k),j(k),l) = DY( ceil(X(i(k),j(k),l)), ceil(Y(i(k),j(k),l)) )...
                             + dy( ceil(X(i(k),j(k),l)), ceil(Y(i(k),j(k),l)) )...
                             .*(Y(i(k),j(k),l)-floor(Y(i(k),j(k),l)));
    end

end


%% Compute the gradient of the flow map, i.e. for each X(i,j) from initial position at T0 to final position after time T

% Shifted matrices for central differencing

% X(i-1,j)
iX = circshift(XX,[1 0 0 0]);

% X(i+1,j)
Xi = circshift(XX,[-1 0 0 0]);

% X(i,j-1)
jX = circshift(XX,[0 1 0 0]);

% X(i,j+1)
Xj = circshift(XX,[0 -1 0 0]);


% Y(i-1,j)
iY = circshift(YY,[1 0 0 0]);

% Y(i+1,j)
Yi = circshift(YY,[-1 0 0 0]);

% Y(i,j-1)
jY = circshift(YY,[0 1 0 0]);

% Y(i,j+1)
Yj = circshift(YY,[0 -1 0 0]);

% Gradient (central differencing)

dXY = nan([size(X) 2]);

dXY(:,:,:,1,1) = (Xi(:,:,:,2)-iX(:,:,:,2))./(Xi(:,:,:,1)-iX(:,:,:,1));

dXY(:,:,:,1,2) = (Xj(:,:,:,2)-jX(:,:,:,2))./(Yj(:,:,:,1)-jY(:,:,:,1));

dXY(:,:,:,2,1) = (Yi(:,:,:,2)-iY(:,:,:,2))./(Xi(:,:,:,1)-iX(:,:,:,1));

dXY(:,:,:,2,2) = (Yj(:,:,:,2)-jY(:,:,:,2))./(Yj(:,:,:,1)-jY(:,:,:,1));

% Compute FTLE field

[i,j] = find( isfinite( squeeze(XX(:,:,1,1)) ) );

FTLE = nan(size(dx,1),size(dx,2),size(XX,3));

for k = 1:length(i)

    for t = 1:size(XX,3)

        if isfinite( squeeze(dXY(i(k),j(k),t,:,:)) )
            FTLE(mi+i(k),ni+j(k),t) = log( max( eig( squeeze(dXY(i(k),j(k),t,:,:))'*squeeze(dXY(i(k),j(k),t,:,:)) ) )^0.5 ) /dt;
        end

    end

end


%% Control plots

%figure

%imagescnan( [], [], flipdim( permute(FTLE(:,:,1),[2 1 3]), 1) )

%colorbar

%figure

%imagescnan( [], [], flipdim( permute(nanmean(FTLE,3),[2 1 3]), 1) )

%colorbar

%toc


%% Save FTLE field

%outfile = input('Output file name ? ');
%outfile = [infile(1:length(infile)-3) '_Anticosti_' num2str(dt) 'd'];
outfile = [infile(1:length(infile)-3) '_' num2str(dt) 'd'];

save(outfile,'infile','dt','lon','lat','time','FTLE')

ret = 0;

end
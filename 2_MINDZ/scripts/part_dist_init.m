% This script create the initial distribution of particles
%
% MapsF 2011

%%
% Define geometry of the grid
m = 197; n = 234;

% Open bathymetry file and get nlayer
load('../run/data/bathymetry.mat'); % bottom (depth), depth (edges of layers), nlayer (# of layers)

bottom = permute(bottom,[2 1]);

nlayer = permute(nlayer,[2 1]);


%% Define a uniform initial distribution of particles

part = zeros(n,m);

part(nlayer>0) = 1;

% Get rid of unwanted wet points

part(1:34,:) = 0; part(n-5:n,:) = 0;

part(:,1:5) = 0; part(:,m-5:m) = 0; 

part(130:end,1:27) = 0; part(190:210,180:end) = 0;


part_total = part';

imagesc(flipdim(part,1))

save('../run/data/part_domain_unif.dat','part_total','-ascii')


%% Define a uniform initial distribution of particles inside the GSL

part = zeros(n,m);

part(nlayer>0) = 1;

part(1:110,:) = 0;

part(110:123,55:95) = 0;

part(:,[1:27 185:end]) = 0;

cabot = 120:185;

for i = 1:length(cabot)
    part(1:floor(cabot(i)*1.07-15),cabot(i)) = 0;
end

part_gulf = part';

figure

imagesc(flipdim(part,1))

save('data/part_gulf_unif.dat','part_gulf','-ascii')


%% Define a uniform initial distribution of particles within the GSL constrained by isobaths

part = zeros(n,m);

part(nlayer>0) = 1;

part(1:110,:) = 0;

part(110:123,55:95) = 0;

part(:,[1:27 185:end]) = 0;

cabot = 120:185;

for i = 1:length(cabot)
    part(1:floor(cabot(i)*1.07-15),cabot(i)) = 0;
end

part(bottom<100) = 0; part(bottom>300) = 0;

part_gulf_100 = part';

figure

imagesc(flipdim(part_gulf_100',1))

%save('data/part_gulf_100_300.dat','part_gulf_100','-ascii')

%% Define a uniform initial distribution of particles within the GSL constrained by isobaths 50 m 

part = zeros(n,m);

part(nlayer>0) = 1;

part(1:110,:) = 0;

part(110:123,55:95) = 0;

part(:,[1:27 185:end]) = 0;

cabot = 120:185;

for i = 1:length(cabot)
    part(1:floor(cabot(i)*1.07-15),cabot(i)) = 0;
end

part(bottom<50) = 0; 

part_gulf_50 = part';

figure

imagesc(flipdim(part_gulf_50',1))

save('../run/data/part_gulf_50.dat','part_gulf_50','-ascii')


%% Define an initial distribution of particles in arbitrary patches within the GSL

part = zeros(n,m);

part(nlayer>0) = 1;

% Compute 2D depth field
[i,j] = find(nlayer>0);

depth2d = zeros(size(nlayer));

for k = 1:length(i)
    depth2d(i(k),j(k)) = depth(nlayer(i(k),j(k))+1);
end

% Define base area for patches

part_patch = zeros(n,m);

part_patch(200:219,60:68) = 1; % Patch #1 = Pointe-des-Monts

part_patch(196:210,108:115) = 1; % Patch #2 = Jacques-Cartier Strait

for i = 1:8
    part_patch((185:195)-i,101+i) = 1; % Patch #3 = Honguedo Strait
end

part_patch(nlayer==0) = 0;

figure

imagesc(flipdim(part_patch,1))


% Generic particles

patch = zeros(size(part_patch));

dumm1 = part_patch;

dumm2 = depth2d-130;

dumm1(dumm2<0) = 0;

% Patches along the shelf at their prefered depth

for j = 60:68

    i = find(dumm1(1:219,j)==1,1,'last');

    patch(i,j) = 1;

end

for i = 200:219

    j = find(dumm1(i,1:68)==1,1,'first');

    patch(i,j) = 1;

end

for j = 101:109

    i = find(dumm1(1:195,j)==1,1,'last');

    patch(i,j) = 1;

end


for j = 108:115

    i = find(dumm1(196:210,j)==1,1,'last');

    patch(195+i,j) = 1;

end


% Check that there will be the same # of particles in each patch (base = 1e3)

sum1 = sum(sum(patch(200:219,60:68))); % Patch #1 = Pointe-des-Monts

toto = zeros(size(patch));

toto(200:219,60:68) = patch(200:219,60:68);

patch(toto==1) = 1e3/sum1;


sum2 = sum(sum(patch(196:210,108:115))); % Patch #2 = Jacques-Cartier Strait

toto = zeros(size(patch));

toto(196:210,108:115) = patch(196:210,108:115);

patch(toto==1) = 1e3/sum2;


sum3 = sum(sum(patch(160:195,101:109))); % Patch #3 = Honguedo Strait

toto = zeros(size(patch));

toto(160:195,101:109) = patch(160:195,101:109);

patch(toto==1) = 1e3/sum3;


figure

imagesc(flipdim(depth2d,1))
 
colorbar

caxis([0 130])

hold on

[i,j] = find(nlayer==0);

plot(j,n-i+1,'w.')

[i,j] = find(patch>0);

plot(j,n-i+1,'r+')

plot([58 58],[234-194 234-186],'r')

plot([95 95],[234-199 234-181],'r')

plot([101 101],[234-208 234-201],'r')

axis([55 125 20 55])


patch = patch';

save('data/part_patches.dat','patch','-ascii')


%% Initial distribution at Old Harry (exploratory well site)

part = zeros(n,m);

part(161,140) = 1;

part = part';

save('data/Old_Harry.dat','part','-ascii')


%% Define an initial distribution of particles based on acoustic hotspots

part = zeros(n,m);

% Open hotspot xpo,ypo file
coord = load('../run/data/xy_hotspots.txt');

ncoord = size(coord,1);

for i = 1:ncoord
    part( ceil(coord(i,2)), ceil(coord(i,1))) = 1;
end

part_hot = part';

figure

imagesc(flipdim(part_hot',1))

save('../run/data/part_hotspots.dat','part_hot','-ascii')


% This routine produces the input_particles file from which initial position
% and timing (starting | ending date) of lagrangian particles are read.
% Based on OPA's grid.
%
% MapsF 2011
%


%% Load data from the GSL grid initialization generated with krill_dist_init.m

domain = input('Which domain to initialize? (total|gulf|patch|hot) ');

if strcmp(domain,'total')
    kini = load('../run/data/part_domain_unif.dat');
elseif strcmp(domain,'gulf')
    kini = load('../run/data/part_gulf_unif.dat');
elseif strcmp(domain,'patch')
    kini = load('../run/data/part_patches.dat');
elseif strcmp(domain,'hot')
    kini = load('../run/data/part_hotspots.dat');total
elseif strcmp(domain,'harry')
    kini = load('../run/data/Old_Harry.dat');total
else
    stop('Problem with domain chosen!')
end

% Only consider selected grid cells

[xpo,ypo] = find(kini>0);       % integer indices

xpo = xpo-0.5;                  % real cell center coordinates

ypo = ypo-0.5;

% Add as much particles as needed inside each cell

raninit = input('Random initial position of particles? (true|false) ');

if raninit

    % Random initial horizontal position

    knum = kini(kini>0)*100; % number of particles per cell

    dummx = []; dummy = [];

    for i = 1:length(knum)

        dummx = [dummx; xpo(i)+rand(floor(knum(i)),1)-0.5]; 

        dummy = [dummy; ypo(i)+rand(floor(knum(i)),1)-0.5];

    end

    xpo = dummx;

    ypo = dummy;

end


%% Construct array in the appropriate format

zpo = ones(length(xpo),1);      % vertical position

days = 30+zeros(length(xpo),1); % duration of advection, in days

year = input('Start year? ');

back = input('Backward simulation? (true|false) ');

if back
    tini = [0 0 0 31 10 year];      % starting date of particle advection

    dini = datenum(year,10,31);     % first date, matlab format
else
    tini = [0 0 0 1 3 year];
    %tini = [0 0 0 1 1 year];

    dini = datenum(year,3,1);
    %dini = datenum(year,1,1);
end

part = [repmat(tini,length(xpo),1) xpo ypo zpo days];


%% Add particles at specific time(s)

for t = 2:2:244 %2:365
    
    if back
        tini = flipdim( datevec(dini-t), 2);
    else
        tini = flipdim( datevec(dini+t), 2);
    end
    
    dummy = [repmat(tini,length(xpo),1) xpo ypo zpo days];

    part = [part; dummy];
    
end


%% Save array in an ascii file in the appropriate format

outfile = input('output file name? ');

fid = fopen(['/workspace/users/OPA/maps/IBMS/Krill/Run/init/' outfile],'w+');

fprintf(fid,'%1i %1i %1i %2i %2i %4i %5.1f %5.1f %5.1f %2i\n',part');

fclose(fid);

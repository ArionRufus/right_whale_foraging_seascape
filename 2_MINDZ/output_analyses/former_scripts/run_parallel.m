% This routine allows the parallelisation of a MATLAB routine.
%
% MapsF 2014

tic

%% Application-specific code

infile = input('Input file list ? ');

fid = fopen(infile,'r');
list = textscan(fid,'%s');
fclose(fid);

ll = length(list{1,1});

dt = input('Time step for computing FTLE (in days, even)? ');

%% Run parallelized routine

core = min(12,ll);

par = parpool(core);

spmd

    for i=drange(1:ll)

        infile = list{1}{i};

        ret = FTLE_comp_back(infile,dt);

    end

end

toc

delete(par)

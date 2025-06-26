% This routine extracts FTLE results used to produce a climatology.
%
% MapsF 2013

%  load data

infirst = input('Enter input filename for the first year of FTLE: ');

year    = input('Enter first and end years for data extraction [y1,y2]: ');

month   = input('Enter the month for data extraction [real]: ');


%% Extract & save data

torigin = datenum(2005,5,1);

dummy = [];

iy = strfind(infirst,num2str(year(1)));

for y = year(1):year(2)
    
    infile = [infirst(1:iy-1) num2str(y) infirst(iy+4:end)];
    
    load(infile)

    t = time>=datenum(y,month,1)-torigin & time<=datenum(y,month,31)-torigin;

    dummy = cat(3,dummy,FTLE(:,:,t));

end

FTLE = dummy;

outfile = [infirst(1:iy-1) num2str(year(1)) '-' num2str(year(2)) '_' num2str(month) infirst(iy+4:end)];

save(outfile,'FTLE','lon','lat','dt')

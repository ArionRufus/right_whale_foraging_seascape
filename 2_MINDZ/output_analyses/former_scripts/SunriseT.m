function [SunR,SunS] = SunriseT(day,lon,lat)
%
% Calcul de l'heure de lever ou coucher du soleil selon la date et la latitude.
% Retourne SunR & SunS en fraction de jour.
% Attention au fuseau horaire!
%
% Inputs:
%  day = Julian day
%  lon = longitude
%  lat = latitude
%
% Outputs:
%  SunR = Sun rise time, in fraction of a day
%  SunS = Sun set time, in fraction of a day
%
% Modif. MapsF 2011


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tzone = -5; %time zone, + pour E

% Four eq below in radian
gamma = 2.*pi./365 .* floor(day);

eqtime = 229.18 .* (  0.000075 + 0.001868.*cos(gamma) ...
                    - 0.032077.*sin(gamma) ...
                    - 0.014615.*cos(2.*gamma) ...
                    - 0.040849.*sin(2.*gamma) );

decl =  0.006918 - 0.399912.*cos(gamma) + 0.070257.*sin(gamma) ...
      - 0.006758.*cos(2.*gamma) + 0.000907.*sin(2.*gamma) ...
      - 0.002697.*cos(3.*gamma) + 0.00148.*sin(3.*gamma);

% !!! ATTENTION: mix of degree & radian !!!
ha = acos( cosd(90.833)./(cosd(lat).*cos(decl)) - tand(lat).*tan(decl) );

% Convert from radian to degree
ha = ha.*180./pi;

SunR = mod( round( 720+4.*(lon-ha)-eqtime )-tzone*60, 1440 )/60;

SunS  = mod( round( 720+4.*(lon+ha)-eqtime )-tzone*60, 1440 )/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resultats

if 0==1

    fid = fopen('sunrise_sunset_2008.txt','w');

    fprintf(fid,'\n Latitude: %5.2f, Longitude: %5.2f\n',lat,lon);

    fprintf(fid,'\nDate          Lever     Midi loc  Coucher \n');

    for ii = 1:length(SunRTime)
    %    for ii = 150:155
        fprintf(fid,'%s   %s  %s  %s\n',datestr(DayofY(ii),1),datestr(SunR(ii),13),datestr(SolNoon(ii),13),datestr(SunS(ii),13));
    end

    fclose(fid);

end

end
function Map=areaweighting(latres,lonres,polarbox,flag)
% Map=areaweighting(latres,lonres,polarbox,flag)
%
% Ben Kravitz (ben.kravitz.work@gmail.com)  12 October 2011
%
% This function creates a map of area weightings of the globe divided
% into grid boxes.  lonres and latres are the resolution of the grid
% in degrees.  (An example is latres=4 and lonres=5 for what is commonly
% known in climate modeling as 4x5 resolution.)
%
% If polarbox==1, the poles have half-boxes, i.e. the pole is the center
% of a grid "box".  If polarbox==0, the poles do not, i.e., the pole is
% the edge of a grid "box".  If one of these settings doesn't match your
% resolution, try the other one.  There's something about the way I
% coded this that makes it work well for some resolutions and poorly for
% others.  I could fix this, but to be honest, I don't feel like
% it right now.  If you want, please feel free to correct the code and
% email it to me.  I'll happily give you credit.
%
% If flag==1, radius depends on latitude. If flag==0, the mean Earth
% radius is used.  Your answer won't change much between these two
% settings, but it's a fun feature that wasn't hard to add.

if polarbox==0;
    latrange=linspace(-90,90,180/latres+1);
    latbottoms=latrange(1:length(latrange)-1);
else
    numlats=180/latres+1;
    latbottoms=(-90)*ones(numlats,1);
    latbottoms(2:numlats)=linspace(-90+latres/2,90-latres/2,180/latres);
end

a=6378.1370; % Earth's semimajor axis (km)
b=6356.7523; % Earth's semiminor axis (km)
rmean=6371.0; % mean radius (km)
if flag==1;
    rnumer=((a^2)*cos(latbottoms)).^2 + ((b^2)*sin(latbottoms)).^2;
    rdenom=(a*cos(latbottoms)).^2 + (b*sin(latbottoms)).^2;
    r=((rnumer)./(rdenom)).^(0.5); % latitude-dependent earth radius
else
    r=rmean*ones(length(latbottoms),1);
end

lonrange=linspace(0,360-lonres,360/lonres);
lat=length(latbottoms);
latusable=zeros(lat+1,1);
latusable(1:lat)=latbottoms;
latusable(lat+1)=90;
areaout=zeros(lat,1);

rmod=rmean*ones(length(r)+1,1);
rmod(1:length(r))=r;
if flag==1;
    rmod(length(r)+1)=a;
end
ruse=zeros(lat,1);
colats=zeros(lat+1,1);
dphi=lonres*pi/180;
ruse(1:lat)=(rmod(1:lat)+rmod(2:lat+1))/2; % mean radius for each lat
colats(1:lat)=pi/2-latusable(1:lat)*pi/180;
areaout(1:lat)=abs((cos(colats(2:lat+1))-cos(colats(1:lat)))*dphi.*(ruse.^2));
biggest=max(areaout);
mapcols=areaout/biggest;
Map=repmat(mapcols,[1 length(lonrange)]);
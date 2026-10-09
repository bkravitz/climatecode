function [Map,Area]=areaweighting(latres,lonres,polarbox,flag)
% AREAWEIGHTING  Map of relative grid-box areas for a regular lat-lon grid
%
% Map=areaweighting(latres,lonres,polarbox,flag)
% [Map,Area]=areaweighting(latres,lonres,polarbox,flag)
%
% Ben Kravitz (ben.kravitz.work@gmail.com)  12 October 2011
%
% This function creates a map of area weightings of the globe divided
% into grid boxes.  lonres and latres are the resolution of the grid
% in degrees.  (An example is latres=4 and lonres=5 for what is commonly
% known in climate modeling as 4x5 resolution.)
%
% The area of a box between latitudes lat1 and lat2 that is dlon wide is
%   r^2 * dlon * |sin(lat2) - sin(lat1)|   (dlon in radians),
% which is exact on a sphere of radius r.
%
% Inputs:
%   latres   - latitude spacing (degrees).  180/latres should be an
%              integer; otherwise it is rounded to the nearest whole
%              number of boxes (with a warning).
%   lonres   - longitude spacing (degrees), same rule for 360/lonres
%   polarbox - [optional, default 0]
%              If polarbox==1, the poles have half-boxes, i.e. the pole is
%              the center of a grid "box" (e.g. CESM finite volume grids,
%              180/latres+1 latitudes, the polar boxes latres/2 tall).
%              If polarbox==0, the pole is the edge of a grid "box"
%              (180/latres latitudes, all latres tall).
%   flag     - [optional, default 0]
%              If flag==1, radius depends on latitude (geocentric radius of
%              the WGS84 ellipsoid, averaged over the box's two latitude
%              edges; an approximation to the true ellipsoidal area).
%              If flag==0, the mean Earth radius (6371 km) is used.  Your
%              answer won't change much between these two settings, but
%              it's a fun feature that wasn't hard to add.
%
% Outputs:
%   Map  - nlat x nlon relative areas.  Rows go from south to north (row 1
%          is the box touching -90), columns from longitude 0 eastward.
%          Normalized so the largest box is 1; Map/sum(Map(:)) gives
%          weights that sum to 1.
%   Area - [optional] the same box areas in km^2 (sums to 4*pi*6371^2 for
%          flag==0)
%
% Bug fixes 8 October 2026 (with Claude): flag==1 took cos/sin of
%   latitudes in degrees as if radians and used the equatorial radius at
%   the North Pole (now the polar radius), so flag==1 weights were wrong.
%   Resolutions that don't divide 180 (or 360) evenly are now rounded to
%   the nearest number of boxes (previously linspace silently truncated
%   for polarbox==0 and polarbox==1 crashed).  polarbox and flag are now
%   optional, and the second output Area was added.

if nargin<3 || isempty(polarbox)
    polarbox=0;
end
if nargin<4 || isempty(flag)
    flag=0;
end
if ~isscalar(latres) || ~isscalar(lonres) || latres<=0 || lonres<=0
    error('areaweighting: latres and lonres must be positive scalars.');
end

% number of boxes (whole numbers)
nlat=round(180/latres);
nlon=round(360/lonres);
if abs(180/latres-nlat)>1e-6 || abs(360/lonres-nlon)>1e-6
    warning('areaweighting:resolution',['180/latres or 360/lonres is not ' ...
        'an integer; using %d latitude spacings and %d longitudes.'],nlat,nlon);
end

% latitude edges of the boxes (degrees), south to north
if polarbox==0
    latedges=linspace(-90,90,nlat+1)';
else
    dlat=180/nlat;
    latedges=[-90 ; linspace(-90+dlat/2,90-dlat/2,nlat)' ; 90];
end

a=6378.1370; % Earth's semimajor axis (km)
b=6356.7523; % Earth's semiminor axis (km)
rmean=6371.0; % mean radius (km)
if flag==1
    rnumer=((a^2)*cosd(latedges)).^2 + ((b^2)*sind(latedges)).^2;
    rdenom=(a*cosd(latedges)).^2 + (b*sind(latedges)).^2;
    redge=sqrt(rnumer./rdenom); % latitude-dependent (geocentric) earth radius
    ruse=(redge(1:end-1)+redge(2:end))/2; % mean radius for each box
else
    ruse=rmean*ones(numel(latedges)-1,1);
end

dphi=(360/nlon)*pi/180; % box width (radians)
Area=abs(diff(sind(latedges)))*dphi.*(ruse.^2); % box area (km^2), one per latitude
Area=repmat(Area,[1 nlon]);
Map=Area/max(Area(:));

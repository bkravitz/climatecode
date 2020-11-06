function mhf = calc_mhf(flux,lat)

% function outvar = calc_mhf(invar,lat)
%
% This function calculates meridional heat flux from a 1-dimensional
% (zonally averaged, time averaged) flux field.  To calculate total
% meridional heat flux, input the net top-of-atmosphere flux
% (rsdt-rsut-rlut).  For oceanic heat flux, input the net surface fluxes
% (rsds-rsus+rlds-rlus).  For atmospheric heat flux, input the
% difference between the two.
%
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 3 November 2020
%
% Input:  flux = radiative flux (W/m2)
%         lat  = latitudes in the flux file
%
% Output:  mhf = meridional heat flux (W)
%                divide by 1e15 to get PW

% Calculating latitude weights
I=size(lat);
if I(1)>I(2);
    lat=transpose(lat);
end
latedges=unique([-90 lat 90]);
if (length(latedges)-1)~=length(flux);
    latedges=unique([-90 (latedges(1)+latedges(2))/2:(latedges(2)-latedges(1)):90 90]);
end

colatedges=pi/2-latedges*pi/180;
dphi=1*pi/180;
areaout=abs(diff(cos(colatedges)))*dphi.*(6371000)^2;
areaout=repmat(areaout,[360 1]);
latwts=sum(areaout,1);

I=size(flux);
if I(1)>I(2);
    flux=transpose(flux);
end

fluxavg=sum(flux.*latwts)./sum(latwts);

mhf = cumsum((flux-fluxavg).*latwts);

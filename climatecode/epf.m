function [vr,wr,srf,epy,epz,epd,vptp100,vpthetap,theta,epd_acc]=epf(t,q,u,v,w,wflag,lats,lons,press)

% EPF  Eliassen-Palm flux and Transformed Eulerian Mean (TEM) circulation
%
% [vr,wr,srf,epy,epz,epd,vptp100,vpthetap,theta,epd_acc] = ...
%     epf(t,q,u,v,w,wflag,lats,lons,press)
%
% Uses the log-pressure, spherical TEM equations of Andrews, Holton &
% Leovy (1987), "Middle Atmosphere Dynamics", eqs. 3.5.1-3.5.5:
%
%   v* = vbar - (1/rho0) d/dz( rho0 v'th'/th_z )
%   w* = wbar + 1/(a cos(phi)) d/dphi( cos(phi) v'th'/th_z )
%   F_phi = rho0 a cos(phi) ( u_z v'th'/th_z - u'v' )
%   F_z   = rho0 a cos(phi) ( [f - 1/(a cos(phi)) d/dphi(ubar cos(phi))] v'th'/th_z - u'w' )
%   div F = 1/(a cos(phi)) d/dphi( cos(phi) F_phi ) + dF_z/dz
%
% where z = -H ln(p/p0) is log-pressure height, rho0 = p/(g H) is the
% matching basic-state density, and w is the log-pressure vertical velocity
% Dz/Dt = -H omega/p.
%
% ---------------------------------------------------------------------
% INPUTS
%
% 3-D arrays, all nlons x nlats x nlevs (for multiple times, loop over time):
%   t     - temperature (K)
%   q     - specific humidity (kg/kg). Only used to get air density when
%           converting w (m/s) to omega (wflag==1). Pass [] for dry air.
%   u     - zonal wind (m/s)
%   v     - meridional wind (m/s)
%   w     - vertical velocity:
%             wflag==1: w is geometric vertical velocity (m/s)
%             wflag==2: w is omega (Pa/s)
%
% 1-D vectors (row or column):
%   lats  - latitudes, degrees or radians, in any order. Values with
%           max(abs(lats)) > pi/2 are treated as degrees.
%   lons  - longitudes (only their number is used, to check array sizes)
%   press - pressure at each level, Pa, in any order. If max(press) <= 1100
%           the values are assumed to be hPa and are converted, with a warning.
%
% ---------------------------------------------------------------------
% OUTPUTS
%
% 2-D zonal-mean fields, all nlats x nlevs:
%   vr      - residual meridional velocity v* (m/s)
%   wr      - residual vertical velocity w* (m/s, log-pressure)
%   srf     - residual mass stream function (kg/s),
%             psi* = 2 pi a cos(phi) * integral from z to model top of rho0 v* dz'
%             so that rho0 v* = -(1/(2 pi a cos(phi))) dpsi*/dz. It is set
%             to zero at the highest model level, so it is only meaningful
%             if that level is near the top of the atmosphere.
%   epy     - meridional EP-flux component F_phi (kg/s^2)
%   epz     - vertical EP-flux component F_z (kg/s^2)
%   epd     - EP-flux divergence div F (kg/m/s^2). Positive means
%             divergence (westerly forcing), negative means convergence.
%   epd_acc - zonal wind forcing by the eddies, div F/(rho0 a cos(phi)),
%             in m/s/day
%
% Other outputs:
%   vptp100  - v'T' at the level nearest 100 hPa (K m/s). This is the full
%              nlons x nlats field, NOT a zonal mean.
%   vpthetap - v'theta' (K m/s), nlons x nlats x nlevs
%   theta    - potential temperature (K), nlons x nlats x nlevs
%
% NaN is returned wherever the result can't be computed: missing input
% data (e.g. levels below ground), the poles (where 1/cos(phi) is
% infinite) and points where the stratification d(theta)/dz is zero.
% EP-flux vectors are returned in physical units, unscaled. For plotting,
% apply your own scaling (e.g. Edmon et al. 1980, or Jucker 2021).
%
% Requires deriv.m (3-point derivative on a non-uniform grid).

% ---------------------------------------------------------------------
% check inputs

if wflag~=1 && wflag~=2
    error('epf: wflag must be 1 (w in m/s) or 2 (w is omega in Pa/s).');
end

nlons=numel(lons);
nlats=numel(lats);
nlevs=numel(press);

expected=[nlons nlats nlevs];
names={'t','u','v','w'};
fields={t,u,v,w};
for k=1:numel(fields)
    if ~isequal(size(fields{k}),expected)
        error('epf: %s must be nlons x nlats x nlevs = %d x %d x %d.',names{k},nlons,nlats,nlevs);
    end
end
if isempty(q)
    q=zeros(expected);
elseif ~isequal(size(q),expected)
    error('epf: q must be nlons x nlats x nlevs, or [].');
end
if nlats<3 || nlevs<3
    error('epf: need at least 3 latitudes and 3 levels to take derivatives.');
end

% latitude as a column (nlats x 1) and pressure as a row (1 x nlevs), so
% that 2-D zonal-mean fields (nlats x nlevs) broadcast via implicit expansion
lat=lats(:);
press=press(:)';

% convert latitude to radians if it is in degrees
if max(abs(lat))>pi/2
    lat=lat*pi/180;
end

% convert pressure to Pa if it looks like hPa
if max(press)<=1100
    warning('epf:pressUnits',['max(press) = %g, so assuming press is in hPa ' ...
        'and converting to Pa. Pass Pa to avoid this warning.'],max(press));
    press=press*100;
end

% ---------------------------------------------------------------------
% constants

rd=287;          % gas constant for dry air (J/kg/K)
rh2o=461;        % gas constant for water vapor (J/kg/K)
cp=1004;         % specific heat of dry air at constant pressure (J/kg/K)
kappa=rd/cp;     % R/cp
g=9.81;          % gravitational acceleration (m/s^2)
a=6.371e6;       % Earth mean radius (m)
p0=1e5;          % reference pressure (Pa)
Omega=7.292e-5;  % Earth's rotation rate (rad/s)
H=rd*240/g;      % log-pressure scale height (about 7 km, for T = 240 K)

% ---------------------------------------------------------------------
% coordinates and basic state

f=2*Omega*sin(lat);    % Coriolis parameter (nlats x 1)
coslat=cos(lat);       % for multiplying (correctly goes to 0 at the poles)

% for dividing: NaN at the poles, where 1/cos(phi) is infinite
coslat_div=coslat;
coslat_div(abs(abs(lat)-pi/2)<1e-6)=NaN;

z=-H*log(press/p0);    % log-pressure height (m), 1 x nlevs
rho0=press/(g*H);      % basic-state density, = (p0/(g H)) exp(-z/H), 1 x nlevs

p3=reshape(press,1,1,nlevs);    % pressure as 1 x 1 x nlevs, for 3-D fields

% ---------------------------------------------------------------------
% thermodynamics and vertical velocity

theta=t.*(p0./p3).^kappa;       % potential temperature (K)

% convert vertical velocity to omega (Pa/s) if needed, using the actual
% moist air density from the ideal gas law (omega = -rho g w)
if wflag==1
    rgas=rd*(1-q)+rh2o*q;       % gas constant for moist air
    dens=p3./(rgas.*t);         % air density (kg/m^3)
    oma=-dens*g.*w;
else
    oma=w;
end

% log-pressure vertical velocity Dz/Dt = -H omega/p (m/s). This is the
% vertical velocity that goes with the z coordinate used for derivatives.
wlp=-H*oma./p3;

% ---------------------------------------------------------------------
% zonal means and deviations from them

% zonal mean, reshaped to nlats x nlevs (no squeeze, so it is safe when a
% dimension has length 1)
zm=@(x) reshape(mean(x,1,'omitnan'),nlats,nlevs);

thetaz=zm(theta);
uz=zm(u);
vz=zm(v);
wz=zm(wlp);

% deviations from the zonal mean (implicit expansion along longitude)
tp=t-mean(t,1,'omitnan');
thetap=theta-mean(theta,1,'omitnan');
up=u-mean(u,1,'omitnan');
vp=v-mean(v,1,'omitnan');
wp=wlp-mean(wlp,1,'omitnan');

vptp=vp.*tp;          % v'T'
vpthetap=vp.*thetap;  % v'theta'

% zonal-mean eddy fluxes (nlats x nlevs)
vpthetapz=zm(vpthetap);   % [v'theta']
upvpz=zm(up.*vp);         % [u'v']
upwpz=zm(up.*wp);         % [u'w']

% ---------------------------------------------------------------------
% TEM residual circulation

dthdz=deriv(thetaz,z,2);    % d(theta)/dz
dudz=deriv(uz,z,2);         % d(ubar)/dz

% v'theta'/theta_z (m^2/s), which appears in every TEM term
psi_e=vpthetapz./dthdz;

vr=vz-deriv(rho0.*psi_e,z,2)./rho0;                  % v* (m/s)
wr=wz+deriv(coslat.*psi_e,lat,1)./(a*coslat_div);    % w* (m/s)

% residual mass stream function: integrate rho0 v* down from the model top.
% Sort levels so z increases, integrate from the top down (flip, then
% cumtrapz), and put the levels back in their original order.
[zs,iz]=sort(z);
integrand=rho0(iz).*vr(:,iz);
integral_top=-flip(cumtrapz(flip(zs),flip(integrand,2),2),2);   % = int_z^top
srf=zeros(nlats,nlevs);
srf(:,iz)=2*pi*a*coslat.*integral_top;

% ---------------------------------------------------------------------
% EP flux

% 1/(a cos(phi)) d/dphi(ubar cos(phi)), the relative vorticity part of
% the absolute vorticity
dudlat=deriv(coslat.*uz,lat,1)./(a*coslat_div);

epy=rho0.*a.*coslat.*(dudz.*psi_e-upvpz);             % F_phi
epz=rho0.*a.*coslat.*((f-dudlat).*psi_e-upwpz);       % F_z

% EP-flux divergence
epd=deriv(coslat.*epy,lat,1)./(a*coslat_div)+deriv(epz,z,2);

% zonal wind forcing by the eddies, converted to m/s/day
epd_acc=epd./(rho0.*a.*coslat_div)*86400;

% ---------------------------------------------------------------------
% flag points that can't be computed (Inf from theta_z = 0, poles, etc.)

vr(~isfinite(vr))=NaN;
wr(~isfinite(wr))=NaN;
srf(~isfinite(srf))=NaN;
epy(~isfinite(epy))=NaN;
epz(~isfinite(epz))=NaN;
epd(~isfinite(epd))=NaN;
epd_acc(~isfinite(epd_acc))=NaN;

% ---------------------------------------------------------------------
% v'T' at the level nearest 100 hPa (full lon x lat field)

[~,k100]=min(abs(press-1e4));
vptp100=vptp(:,:,k100);

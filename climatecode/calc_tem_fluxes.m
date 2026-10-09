function [vstar,wstar] = calc_tem_fluxes(lat,lev,theta_bar,v_bar,w_bar,vth_bar,ccmi_flag,omega_flag,varargin)

% CALC_TEM_FLUXES  TEM residual velocities v* and w* from zonal means
%
% [vstar,wstar] = calc_tem_fluxes(lat,lev,theta_bar,v_bar,w_bar,vth_bar,ccmi_flag,omega_flag,[wv])
%
% This script calculates the residual velocities v-star and w-star.
% Script originally written by Rolando Garcia
% Updated and Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% and Daniele Visioni (dv224@cornell.edu)
% Last updated 13 October 2020
%
% Method: the Transformed Eulerian Mean equations of Andrews, Holton &
% Leovy (1987), "Middle Atmosphere Dynamics", eqs. 3.5.1-3.5.2:
%
%   v* = vbar - (1/rho) d/dz( rho v'th'/th_z )
%   w* = wbar + 1/(a cos(phi)) d/dphi( cos(phi) v'th'/th_z )
%
% Derivatives use 3-point Lagrangian differences (second order on
% non-uniform grids, one-sided at the end points), the same scheme as
% deriv.m.  Latitudes and levels may be in either order.  With
% ccmi_flag = 1 this gives the same v*, w* as epf.m (for the same H).
%
%
%%%%% Inputs %%%%%
%
% lat = latitude (degrees), centers of grid boxes, monotonic (ny values)
% lev = pressure (mb = hPa), monotonic (nz values)
%
% Zonally averaged quantities (ny x nz)
%   theta_bar = potential temperature (K), reference pressure 1000 mb
%   v_bar = meridional velocity (m/s)
%   w_bar = vertical velocity (m/s), or omega (Pa/s) if omega_flag = 1
%   vth_bar = meridional eddy flux of potential temperature [v'theta'] (K m/s)
%
% ccmi_flag  = 1 if using the CCMI standard method of calculating height
%                (log-pressure height z = H ln(1013.25/lev) with H = 7 km
%                and basic-state density rho = p/(g H), as in AHL and
%                epf.m).  w is then the log-pressure vertical velocity.
%            = 0 if computing geometric height directly from temperature
%                (the hypsometric equation, using virtual temperature) and
%                using the actual density p/(R_d T_v).  Horizontal
%                derivatives are still taken along pressure levels.
% omega_flag = 1 if w_bar is omega (Pa/s); it is converted to
%                w = -H omega/p (log-pressure, ccmi_flag = 1) or
%                w = -omega/(rho g) (geometric, ccmi_flag = 0)
%            = 0 if w_bar is already w (m/s)
%
% [optional] wv = zonally averaged water vapor mass mixing ratio (kg/kg),
%                 ny x nz.  Default zero.  Only used for ccmi_flag = 0.
%
%
%%%%% Outputs %%%%%
%
% vstar = residual meridional velocity (m/s), ny x nz
% wstar = residual vertical velocity (m/s), ny x nz
%
% NaN is returned where the result can't be computed (exactly at the
% poles, where d(theta)/dz = 0, or NaN inputs).
%
%
% Note that if you are making calculations for the stratosphere (where
% the TEM velocities make sense to calculate), water vapor ranges
% between about 1 and 7 ppmv, which translates to mass mixing ratios of
% about 6e-7 to 4e-6 kg/kg.  Excluding water vapor therefore changes the
% virtual temperature (and hence density and height) by less than 0.001%.
% As such, although we allow for the option of water vapor mass
% mixing ratios to be included, it is not necessary.  The mass mixing
% ratio of water vapor is set to zero by default.
%
% Mastenbrook, H. J. (1974), Water-vapor measurements in the lower
% stratosphere, Canadian Journal of Chemistry, 52, 1527-1531,
% doi:10.1139/v74-224.
%
% Bug fixes 8 October 2026 (with Claude):
%   - w* eddy term was dimensionally wrong (d/dy of v'th'/th_z divided by
%     a band width, with no cos(phi) inside the derivative) and the first
%     and last latitudes were overwritten with an ad hoc pole formula; now
%     AHL eq. 3.5.2.  This changes w* everywhere.
%   - omega -> w for ccmi_flag = 1 used the actual density (geometric w);
%     it is now the log-pressure w = -H omega/p (AHL section 3.1), which
%     is what goes with the log-pressure z.  Changes w* by up to ~20%.
%   - ccmi_flag = 0: geometric heights had the wrong sign (or were
%     computed downward from the first level), were shifted half a level,
%     used an invented level at lev(end)/10, and the array size was hard
%     coded to 192 latitudes (crashed otherwise).  Heights are now the
%     hypsometric heights of the pressure levels themselves.
%   - the near-surface d(theta)/dz floor indexed levels by latitude
%     indices (crashed whenever lev(1) was near the surface); removed.
%   - T from theta used p0 = 1013.25 mb instead of 1000 mb; the moist
%     specific heat used liquid water (4184) instead of vapor; virtual
%     temperature was combined with the moist gas constant (moisture
%     counted twice).  These only affect ccmi_flag = 0.
%   - derivatives were interpolated back to the levels linearly in
%     pressure, giving NaN at the top/bottom levels and the end latitudes;
%     now 3-point derivatives directly on the levels (no NaN edges).
%   - (a+z) replaced by a (shallow atmosphere, as in AHL); debug print of
%     mean(z) removed; lat and lev may now be rows or columns.

% checking inputs
if nargin<8
    error('calc_tem_fluxes: need at least 8 inputs.');
end
lat=lat(:);  % column (ny x 1)
lev=lev(:)'; % row (1 x nz), so ny x nz fields broadcast
ny=numel(lat);
nz=numel(lev);
if ny<3 || nz<3
    error('calc_tem_fluxes: need at least 3 latitudes and 3 levels.');
end
if ~(all(diff(lat)>0) || all(diff(lat)<0)) || ~(all(diff(lev)>0) || all(diff(lev)<0))
    error('calc_tem_fluxes: lat and lev must be strictly monotonic.');
end
names={'theta_bar','v_bar','w_bar','vth_bar'};
fields={theta_bar,v_bar,w_bar,vth_bar};
for k=1:numel(fields)
    if ~isequal(size(fields{k}),[ny nz])
        error('calc_tem_fluxes: %s must be ny x nz = %d x %d.',names{k},ny,nz);
    end
end
if ~isscalar(ccmi_flag) || ~any(ccmi_flag==[0 1]) || ~isscalar(omega_flag) || ~any(omega_flag==[0 1])
    error('calc_tem_fluxes: ccmi_flag and omega_flag must be 0 or 1.');
end

% filling the water vapor array
% only necessary for ccmi_flag = 0
if isempty(varargin) || isempty(varargin{1}) % water vapor was not input into the script
    wv=zeros(ny,nz);
else
    wv=varargin{1};
    if ~isequal(size(wv),[ny nz])
        error('calc_tem_fluxes: wv must be ny x nz = %d x %d.',ny,nz);
    end
end

% Defining constants
H =  7000;     % scale height (m)
ae = 6371000;  % Earth radius (m)
g =  9.80665;  % acceleration due to gravity (m/s2)
ps = 1013.25;  % sea level pressure (mb), reference for heights
p0 = 1000;     % reference pressure for potential temperature (mb)
Rd = 287.0;    % gas constant for dry air (J/kg/K)
Rv = 461.52;   % gas constant for water vapor (J/kg/K)
cpd = 1004;    % specific heat of dry air (J/kg/K)
p = 100*lev;   % pressure (Pa), 1 x nz

% Defining the height of each level and the density
if ccmi_flag == 1 % use scale height (log-pressure coordinates)
    z = repmat(H*log(ps./lev),[ny 1]);
    rho = repmat(p/(g*H),[ny 1]); % basic-state density (kg/m3)
    if omega_flag == 1
        w_bar = -H*w_bar./p; % log-pressure vertical velocity (m/s)
    end
else % actually computing the height of the levels
    T = theta_bar.*(lev/p0).^(Rd/cpd); % temperature from potential temperature
    Tv = T.*(1+wv*Rv/Rd)./(1+wv);      % virtual temperature
    rho = p./(Rd*Tv);                  % air density (kg/m3)
    if omega_flag == 1
        w_bar = -w_bar./(rho*g);
    end
    % hypsometric equation, integrating upward from the highest pressure
    % (assumed isothermal between ps and that level)
    [~,ord] = sort(lev,'descend');
    lnp = log(p(ord));
    Tvs = Tv(:,ord);
    zs = zeros(ny,nz);
    zs(:,1) = (Rd*Tvs(:,1)/g)*(log(100*ps)-lnp(1));
    for k=2:nz
        zs(:,k) = zs(:,k-1)+(Rd*(Tvs(:,k-1)+Tvs(:,k))/(2*g))*(lnp(k-1)-lnp(k));
    end
    z = zeros(ny,nz);
    z(:,ord) = zs; % back to the input level order
end

% v'th'/th_z (m2/s), which appears in both TEM terms
d_thbar_dz = ddx(theta_bar,z,2);
X = vth_bar./d_thbar_dz;

% computing TEM velocities (AHL eqs. 3.5.1-3.5.2)
phi = lat*pi/180;
coslat = cos(phi);
coslat_div = coslat;                         % for dividing:
coslat_div(abs(abs(phi)-pi/2)<1e-6) = NaN;   % NaN exactly at the poles
vstar = v_bar - ddx(rho.*X,z,2)./rho;
wstar = w_bar + ddx(coslat.*X,phi,1)./(ae*coslat_div);

% flag points that can't be computed (e.g. Inf from th_z = 0)
vstar(~isfinite(vstar)) = NaN;
wstar(~isfinite(wstar)) = NaN;

end


function r=ddx(y,x,dim)
% dy/dx along dimension dim (1 or 2) of the 2-D array y, using 3-point
% Lagrangian differences (as in deriv.m).  x is a vector or the size of y.
if dim==2
    y=y.';
    x=x.';
end
n=size(y,1);
if isvector(x)
    x=repmat(x(:),[1 size(y,2)]);
end
a=2:n-1;
x01=x(a-1,:)-x(a,:);
x02=x(a-1,:)-x(a+1,:);
x12=x(a,:)-x(a+1,:);
r=zeros(size(y));
r(a,:)=y(a-1,:).*x12./x01./x02 + y(a,:).*(1./x12-1./x01) - y(a+1,:).*x01./x02./x12;
r(1,:)=y(1,:).*(x01(1,:)+x02(1,:))./(x01(1,:).*x02(1,:)) - y(2,:).*x02(1,:)./(x01(1,:).*x12(1,:)) + y(3,:).*x01(1,:)./(x02(1,:).*x12(1,:));
r(n,:)=-y(n-2,:).*x12(end,:)./(x01(end,:).*x02(end,:)) + y(n-1,:).*x02(end,:)./(x01(end,:).*x12(end,:)) - y(n,:).*(x02(end,:)+x12(end,:))./(x02(end,:).*x12(end,:));
if dim==2
    r=r.';
end
end

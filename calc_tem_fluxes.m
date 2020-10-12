function [vstar,wstar] = calc_tem_fluxes(lat,lev,theta_bar,v_bar,w_bar,vth_bar,ccmi_flag,omega_flag,varargin)

% [vstar,wstar] = calc_tem_fluxes(lat,lev,theta_bar,v_bar,w_bar,vth_bar,ccmi_flag,[wv])
%
% This script calculates the residual velocities v-star and w-star.
% Script originally written by Rolando Garcia
% Updated and Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 8 October 2020
%
%
%%%%% Inputs %%%%%
%
% lat = latitude (degrees), centers of grid boxes
% lev = altitude (mb)
%
% Zonally averaged quantities (ny x nz)
%   theta_bar = potential temperature (K)
%   v_bar = meridional velocity (m/s)
%   w_bar = vertical velocity (m/s)
%   vth_bar = meridional eddy flux of potential temperature (K m/s)
%
% ccmi_flag = 1 if using the CCMI standard method of calculating height
%               (using a scale height)
%           = 0 if computing directly from temperature (the hypsometric
%               equation) and, optionally, water vapor mass mixing ratio
% omega_flag = 1 if using omega (Pa/s) instead of w (m/s)
%            = 0 if using w (m/s)
%
% [optional] wv = zonally averaged water vapor mass mixing ratio (kg/kg)
%
%
%%%%% Outputs %%%%%
%
% vstar = residual meridional velocity (m/s)
% wstar = residual vertical velocity (m/s)
%
%
% Note that if you are making calculations for the stratosphere (where
% the TEM velocities make sense to calculate), water vapor ranges
% between 1 and 4 ppm, which translates to mass mixing ratios of
% 8e-4 to 3e-3.  Even at the maximum water vapor concentration in the
% stratosphere, excluding water vapor amounts to maximum errors of 0.1%.
% As such, although we allow for the option of water vapor mass
% mixing ratios to be included, it is not necessary.  The mass mixing
% ratio of water vapor is set to zero by default.
%
% Mastenbrook, H. J. (1974), Water-vapor measurements in the lower
% stratosphere, Canadian Journal of Chemistry, 52, 1527-1531,
% doi:10.1139/v74-224.

% Defining constants
H =  7000;    % scale height (m)
ae = 6371000; % Earth radius (m)
g =  9.81;    % acceleration due to gravity (m/s2)
p0 = 1013.25; % sea level pressure (mb)

% array lengths and other commonly used things
ny = length(lat);
nz = length(lev);

% filling the water vapor array
% only necessary for ccmi_flag = 0
if isempty(varargin)==1; % water vapor was not input into the script
    wv=zeros(ny,nz);
else
    wv=varargin{1};
end

% Defining the height of the middle of each level
if ccmi_flag == 1; % use scale height
    z = repmat(H*log(p0./lev),[1 ny])'; 
    z_m = repmat(H*log(p0./lev),[1 ny-1])'; 
    rho = repmat(100*lev/(g*H),[1 ny])';
    R = 287.0*(1-wv)+461.52*wv; % ideal gas constant (J/kg/K)
    cp = 1004*(1-wv)+4184*wv; % specific heat (J/kg/K)
    g = 9.80665;
    T = theta_bar.*repmat(lev/p0,[1 ny])'.^(R./cp); % temperature from potential temperature
    rho_T = 100*repmat(lev,[1 ny])'./(R.*T);
    if omega_flag == 1
        w_bar = - w_bar./(rho_T*g);
    end
else % actually computing the height of the levels
    R = 287.0*(1-wv)+461.52*wv; % ideal gas constant (J/kg/K)
    cp = 1004*(1-wv)+4184*wv; % specific heat (J/kg/K)
    g = 9.80665;
    T = theta_bar.*repmat(lev/p0,[1 ny])'.^(R./cp); % temperature from potential temperature
    rho = 100*repmat(lev,[1 ny])'./(R.*T);
    if omega_flag == 1
        w_bar = - w_bar./(rho*g);
    end
    Tv = T.*((wv+0.622)./(0.622*(1+wv))); % virtual temperature
    levplus=[lev' lev(end)/10]';
    levfrac=log(levplus(1:end-1)./levplus(2:end));
    h = ((R.*Tv)/g).*repmat(levfrac,[1 ny])'; % level thickness by hypsometric equation
    zedge = zeros(192,length(lev)+1);
    for k=2:length(lev)+1
        zedge(:,k)=zedge(:,k-1)+h(:,k-1);
    end
    z=(zedge(:,1:end-1)+zedge(:,2:end))/2;
    z_m = z(1:size(z,1)-1,:);
end

% determining y-distance of centers of latitude bands from equator
% and edges (needed for interpolation after taking derivatives)
colatcenters=pi/2-lat*pi/180;
y=repmat(cos(colatcenters),[1 nz]).*(ae+z);
latedges=(lat(1:end-1)+lat(2:end))/2; % boundaries between grid boxes
colatedges=pi/2-latedges*pi/180;
y_minus_one=repmat(cos(colatedges),[1 nz]).*(ae+z_m);
latedges_plus=unique([-90 ; latedges ; 90]);
colatedges_plus=pi/2-latedges_plus*pi/180;
area=repmat(abs(diff(cos(colatedges_plus))),[1 nz]).*(ae+z);

% computing level arrays so we can interpolate back after taking derivatives
lev_minus_one = (lev(1:end-1)+lev(2:end))/2;

% Computing vertical derivative of theta (K/m)
% Then interpolating back up to native vertical grid (differencing reduces length)
d_thbar_dz=diff(theta_bar,1,2)./diff(z,1,2);
I=find(z(:,1)<2000);
d_thbar_dz(:,I,:)=max(d_thbar_dz(:,I,:),0.1*ones(size(d_thbar_dz(:,I,:))));
d_thbar_dz_plus = zeros(ny,nz);
for k=1:ny
    d_thbar_dz_plus(k,:)=interp1(lev_minus_one,d_thbar_dz(k,:),lev);
end

% Computing quantity d1z (kg/m2/s)
d1z = diff(rho.*vth_bar./d_thbar_dz_plus,1,2)./diff(z,1,2);
d1z_plus = zeros(ny,nz);
for k=1:ny
    d1z_plus(k,:)=interp1(lev_minus_one,d1z(k,:),lev);
end

% Computing horizontal derivatives
d1y = diff(vth_bar./d_thbar_dz_plus,1,1)./diff(y,1,1);
d1y_plus = zeros(ny,nz);
for k=1:nz;
    d1y_plus(:,k)=interp1(y_minus_one(:,k),d1y(:,k),y(:,k));
end
d2y = diff(d1y_plus,1,1)./diff(y,1,1);
d2y_plus = zeros(ny,nz);
for k=1:nz;
    d2y_plus(:,k)=interp1(y_minus_one(:,k),d2y(:,k),y(:,k));
end

% computing TEM velocities
vstar = v_bar - d1z_plus./rho;
wstar = w_bar + d1y_plus./area;
wstar(1,:) = w_bar(1,:) + d2y_plus(1,:)./(ae+z(1,:));
wstar(end,:) = w_bar(end,:) - d2y_plus(end,:)./(ae+z(end,:));

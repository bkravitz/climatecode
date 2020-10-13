function WBGT=calcwbgt(Ta,RH,S,p,V,tau,zday,lat)

% function WBGT=calcWBGT(Ta,RH,S,p,V,tau,zday,lat)
% Calculates Wet Bulb Globe Temperature
% Written by Ben Kravitz (ben.kravitz@pnnl.gov or ben.kravitz.work@gmail.com)
% Last updated 15 November 2017
%
% Based on formulas provided by
% Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp, and R. Sharp (2008),
% Modeling the wet bulb globe temperature using standard meteorological
% measurements, Journal of Occupational and Environmental Hygiene, 5, 645-655,
% doi:10.1080/15459620802310770.
%
% Inputs
% Ta = surface air temperature (degrees C)
% RH = relative humidity (percent)
% S = total solar irradiance (W/m2)
% p = pressure (Pa)
% V = wind speed (m/s)
% tau = hour of day (0 is midnight, 12 is noon, etc.)
% zday = Julian day (1 through 365 - no leap days)
% lat = latitude (in degrees)
%
% This script is fully vectorized
% To take advantage of this feature, all inputs must be of the same length
%
% Outputs
% WBGT = wet bulb globe temperature (degrees C)
%
% Parametric sensitivity (+/- 10 pct) was performed for inputs
% Ta=40, RH=20, S=1000, p=100000, V=10, tau=12, zday=182, lat=44

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S0=1367; % Solar constant (W/m2)
sigma=5.67e-8; % Stefan Boltzmann constant (W/m2/K4)
Lv=2.453e6; % Latent heat of vaporization (J/kg)
M_H2O=18; % molecular weight of water vapor (g/mol)
M_air=28.97; % molecular weight of air (g/mol)
Rv=461.51; % ideal gas constant for water vapor (J/kg/K)
Rd=287.0; % ideal gas constant for dry air (J/kg/K)
cp_air=1004; % specific heat capacity of air (J/kg/K)
cp_water=4184; % specific heat capacity of water (J/kg/K)
Df=0.282e-4; % Diffusivity of water in air (m2/s) % 10 pct param variance --> 1 pct answer variance
k=4.358e-3; % Thermal conductivity of air (W/m/K) % 10 pct param variance --> 0.4 pct answer variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_w=0.95; % wick emissivity % 10 pct param variance --> 0.02 pct answer variance
eps_g=0.95; % globe emissivity % 10 pct param variance --> 0.3 pct answer variance
alpha_w=0.4; % albedo of the wick % 10 pct param variance --> 0.6 pct answer variance
alpha_g=0.05; % globe albedo % 10 pct param variance --> 0.1 pct answer variance
alpha_sfc=0.45; % surface albedo % 10 pct param variance --> 0.8 pct answer variance
D=0.007; % WBGT wick diameter (m) % 10 pct param variance --> 0.06 pct answer variance
lproj=D/4/0.0254; % D/4L where L = wick length % 10 pct param variance --> 0.06 pct answer variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations below here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat=lat*pi/180;
Ta=Ta+273.15;
I=find(p<1200);p(I)=p(I)*100; % fixing pressures in case they are accidentally specified in mb

q=spechumid(Ta,RH,p); % specific humidity
q(q<0)=0;
cp=cp_air*(1-q)+cp_water*q; % specific heat capacity
rho=p./(Rv*q+Rd*(1-q))./Ta; % air density
ea=(RH/100)*611*(Lv/Rv).*(1/273.15-1./Ta); % partial pressure of water vapor
ea(ea<0)=0;
eps_a=0.575*(ea/100).^0.143; % emissivity of the atmosphere (empirical formula requires ea in mb)

% Calculating fdir (fraction of S that is direct)
[theta,rsdist]=calcsza(tau,zday,lat); % calculating solar zenith angle
fdir=zeros(size(tau));
I=find(theta<=89.5*pi/180);
Smax=S0*cos(theta)./(rsdist.^2);
Sstar=S./Smax;
fdir2=exp(3-1.34*Sstar-1.65./Sstar);
fdir(I)=fdir2(I);

mu=calcviscosity(18.27e-6,291.15,Ta,120);

Pr=cp.*mu/k; % Prandtl number
Sc=mu./rho/Df; % Schmidt number
Re=rho.*V*D./mu; % Reynolds number
h=(k/D)*0.281*Re.^(1-0.4).*Pr.^(1-0.56); % convective heat transfer coefficient
Nu=2.0+0.6*Re.^0.5.*Pr.^(1/3); % Nusselt number, empirically derived
h2=k/D*Nu; % convective heat transfer coefficient derived from Nu

% Initial guesses
Tw=Ta-(100-RH)/5;
Tg=Ta;

% Calculating Tw and Tg iteratively
dT=30*ones(length(tau),2);
flags=dT/0.02;
while sum(flags(:))>=numel(flags);
    I=find(flags(:,1)>=1);
    J=find(flags(:,2)>=1);
    
    FA=sigma*eps_w*(0.5*(1+eps_a(I)).*Ta(I).^4-Tw(I).^4)+(1-alpha_w)*S(I).*((1-fdir(I))*(1+lproj)+fdir(I).*(tan(theta(I))/pi+lproj)+alpha_sfc);
    ew=611*(Lv/Rv)*(1/273.16-1./Tw(I));
    Twnew=Ta(I)-(Lv./cp(I)).*(M_H2O/M_air).*(Pr(I)./Sc(I)).^0.56.*(ew-ea(I))./(p(I)-ew)+FA./h(I);
    Twnew=0.1*Twnew+0.9*Tw(I);
    dT(I,1)=abs(Tw(I)-Twnew);
    Tw(I)=Twnew;
    
    Tgnew=(0.5*(1+eps_a(J)).*Ta(J).^4-h2(J)/eps_g/sigma.*(Tg(J)-Ta(J))+S(J)/2/eps_g/sigma*(1-alpha_g).*(1+(1/2./cos(theta(J))-1).*fdir(J)+alpha_sfc)).^0.25;
    Tgnew=0.1*Tgnew+0.9*Tg(J);
    dT(J,2)=abs(Tg(J)-Tgnew);
    Tg(J)=Tgnew;
    
    flags=dT/0.02;
end

% Putting everything together
WBGT = 0.7*Tw + 0.2*Tg + 0.1*Ta - 273.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta,rsdist]=calcsza(tau,zday,lat)
% Calculating solar zenith angle
% zenit=sin(lat)*sin(delta)+cos(lat)*cos(delta)*cos(H)
% delta=solar declination angle
% H = hour angle
	
eps=23.5*pi/180; % obliquity (radians)
ecc=0.0178; % eccentricity
sol=355; % Julian day of winter solstice
per=4; % Julian day of perihelion
daypyr=365; % number of days per year
rotper=24; % number of hours per day

% Solar Declination Angle
zday=mod(zday,daypyr);
dayper=2*pi/daypyr;
n=3*daypyr/4-sol;
delta=eps*sin((n+zday)*dayper);

% Hour Angle
tofday=mod(tau,rotper);
rot=tofday*2*pi/rotper;
H=rot+pi;

% Corrections to declination angle based on orbital properties
aphel=per+round(daypyr/2);
zdist=(zday-aphel)/daypyr;
rsdist=(1+ecc*cos(2*pi*zdist)).^2;

% Calculating theta
theta=sin(lat).*sin(delta)./rsdist+cos(lat).*cos(delta)./rsdist.*cos(H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q=spechumid(T,RH,p)
% Calculating specific humidity (q)
Lv=2.453e6; % Latent heat of vaporization (J/kg)
Rv=461.5; % ideal gas constant for water vapor (J/kg/K)
Rd=287.0; % ideal gas constant for dry air (J/kg/K)
es=611*(Lv/Rv)*(1/273.16-1./T);
ws=es*Rd/Rv./(p-es);
w=RH/100.*ws;
q=w./(1+w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mu=calcviscosity(mu0,T0,T,C)
% C is Sutherlands constant for the gaseous material
mu=mu0*(T0+C)./(T+C).*(T/T0).^1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

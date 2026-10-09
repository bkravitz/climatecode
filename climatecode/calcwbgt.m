function WBGT=calcwbgt(Ta,RH,S,p,V,tau,zday,lat)

% function WBGT=calcwbgt(Ta,RH,S,p,V,tau,zday,lat)
% Calculates Wet Bulb Globe Temperature
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 15 November 2017
% Corrected 9 August 2026:
%   (1) saturation/partial vapor pressures now use the exponential
%       Clausius-Clapeyron form (the exp() was missing in three places,
%       which understated e_s several-fold at warm temperatures);
%   (2) calcsza now returns the zenith ANGLE via acos (it previously
%       returned the cosine, which downstream cos()/tan() calls then
%       treated as an angle), and the Sun-Earth distance factor no
%       longer contaminates the geometry (it only scales irradiance);
%   (3) added Liljegren's solar consistency cap (calc_solar_parameters):
%       S is rescaled to at most 0.85 of the top-of-atmosphere
%       irradiance at the current sun angle and zeroed below ~0.5 deg
%       solar elevation. Without the cap, time-mean solar inputs near
%       sunrise/sunset blow up the 1/cos(theta) geometry factors by
%       tens of degrees; fdir is now also clamped to [0, 0.9] as in
%       the reference model.
% Bug fixes 8 October 2026 (with Claude):
%   (1) thermal conductivity of air was 4.358e-3 W/m/K (~6x too small;
%       air is ~0.025), which made Pr ~4 instead of ~0.7 and roughly
%       tripled both the evaporative and radiative terms of Tw.  Now
%       k = (cp + 1.25*Rd)*mu, as in Liljegren's code.  Largest change:
%       natural wet bulb (and WBGT) is several degrees warmer, especially
%       at low RH;
%   (2) the globe now uses the 2-inch (0.0508 m) globe diameter of
%       Liljegren et al. instead of the 7 mm wick diameter;
%   (3) the iteration stopped when the MEAN change was below 0.02 K, not
%       every point's, and stopped after one pass for everything if any
%       input was NaN.  Each point now iterates until its own change is
%       < 0.02 K (max 500 iterations, else NaN with a warning); NaN inputs
%       give NaN output;
%   (4) TOA irradiance was multiplied by (r/a)^2 instead of divided by it
%       (the 0.85 cap was ~7% too high near aphelion, too low near
%       perihelion);
%   (5) wind speed has a floor of 0.13 m/s (Liljegren's minimum); V=0
%       gave division by zero;
%   (6) solar geometry now comes from sza.m (Kepler orbit, exact
%       declination, present-day orbit; cos(zenith) = 0 at night, which
%       is excluded from the tan and 1/cos terms below); see sza.m;
%   (7) inputs may now be any mix of scalars and same-size arrays
%       (previously only same-length vectors worked).
% Bug fix 9 October 2026 (with Claude): a missing (NaN) solar irradiance
%   S gave a finite WBGT (MATLAB's min() ignores NaN, so it was replaced by
%   the 0.85 cap in daylight and 0 at night); it now gives NaN, as stated.
%
% Based on formulas provided by
% Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp, and R. Sharp (2008),
% Modeling the wet bulb globe temperature using standard meteorological
% measurements, Journal of Occupational and Environmental Hygiene, 5, 645-655,
% doi:10.1080/15459620802310770.
%
% Method
% Natural wet bulb temperature Tw (Liljegren eq. 6, cylindrical wick) and
% globe temperature Tg (Liljegren eq. 8, 2-inch black globe) are each
% found by under-relaxed fixed-point iteration (Tnew = 0.9*Told +
% 0.1*Tcalc) until the change is < 0.02 K, as in Liljegren's code.
% Air properties are evaluated at Ta (Liljegren uses the film temperature
% for some); the surface temperature is assumed equal to Ta.  Vapor
% pressures use Clausius-Clapeyron over liquid water with constant Lv.
% Outdoor WBGT = 0.7*Tw + 0.2*Tg + 0.1*Ta.
% Requires sza.m (same folder) for the solar zenith angle.
%
% Inputs (scalars or arrays of the same size; scalars are expanded)
% Ta = surface air temperature (degrees C)
% RH = relative humidity (percent)
% S = total (global horizontal) solar irradiance (W/m2)
% p = pressure (Pa; values < 1200 are assumed to be hPa and converted)
% V = wind speed (m/s), ideally at 2 m; floored at 0.13 m/s
% tau = hour of day (local solar time; 0 is midnight, 12 is noon, etc.)
% zday = Julian day (1 through 365 - no leap days)
% lat = latitude (in degrees)
%
% Outputs
% WBGT = wet bulb globe temperature (degrees C), same size as the inputs
%
% Parametric sensitivity (+/- 10 pct) was performed for inputs
% Ta=40, RH=20, S=1000, p=100000, V=10, tau=12, zday=182, lat=44
% (before the 8 October 2026 fixes; the percentages below may differ now)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_w=0.95; % wick emissivity % 10 pct param variance --> 0.02 pct answer variance
eps_g=0.95; % globe emissivity % 10 pct param variance --> 0.3 pct answer variance
alpha_w=0.4; % albedo of the wick % 10 pct param variance --> 0.6 pct answer variance
alpha_g=0.05; % globe albedo % 10 pct param variance --> 0.1 pct answer variance
alpha_sfc=0.45; % surface albedo % 10 pct param variance --> 0.8 pct answer variance
D=0.007; % WBGT wick diameter (m) % 10 pct param variance --> 0.06 pct answer variance
lproj=D/4/0.0254; % D/4L where L = wick length (0.0254 m) % 10 pct param variance --> 0.06 pct answer variance
Dg=0.0508; % globe diameter (m), 2-inch globe as in Liljegren et al. (2008)
Vmin=0.13; % minimum wind speed (m/s), as in Liljegren's code
tol=0.02; % convergence criterion (K)
maxiter=500; % maximum number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations below here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Expand scalars to a common size and work with column vectors
try
    z=0*(Ta+RH+S+p+V+tau+zday+lat); % zeros of the common size
catch
    error('calcwbgt:size','Inputs must be scalars or arrays of the same size.');
end
sz=size(z);
z(isnan(z))=0;
Ta=Ta(:)+z(:); RH=RH(:)+z(:); S=S(:)+z(:); p=p(:)+z(:);
V=V(:)+z(:); tau=tau(:)+z(:); zday=zday(:)+z(:); lat=lat(:)+z(:);

lat=lat*pi/180;
Ta=Ta+273.15;
I=p<1200; p(I)=p(I)*100; % fixing pressures in case they are accidentally specified in mb
V(V<Vmin)=Vmin; % Liljegren's minimum wind speed (also avoids division by zero)
S(S<0)=0; % (comparisons, unlike max(), keep NaN as NaN)

q=spechumid(Ta,RH,p); % specific humidity
q(q<0)=0;
cp=cp_air*(1-q)+cp_water*q; % specific heat capacity
rho=p./(Rv*q+Rd*(1-q))./Ta; % air density
ea=(RH/100)*611.*exp((Lv/Rv)*(1/273.15-1./Ta)); % partial pressure of water vapor
ea(ea<0)=0;
eps_a=0.575*(ea/100).^0.143; % emissivity of the atmosphere (empirical formula requires ea in mb)

% Solar consistency (Liljegren calc_solar_parameters): rescale S to at
% most 0.85 of the top-of-atmosphere irradiance at the current sun
% angle, zero it below ~0.5 deg solar elevation, and compute the
% direct-beam fraction from the rescaled ratio.
[coszen,fluxfac]=sza(tau,zday,lat); % cos(zenith) (night = 0) and (a/r)^2
theta=acos(min(coszen,1)); % solar zenith ANGLE (radians)
Smax=S0*coszen.*fluxfac; % top-of-atmosphere irradiance (W/m2)
Smax(theta>89.5*pi/180)=0;
Snan=isnan(S); % remember missing irradiance before capping
Sstar=zeros(size(Ta));
I=find(Smax>0);
Sstar(I)=min(S(I)./Smax(I),0.85);
S=Sstar.*Smax; % capped (and night-zeroed) solar irradiance
S(isnan(Sstar) | Snan)=NaN; % missing irradiance stays missing (min() above ignores NaN)
fdir=zeros(size(Ta)); % direct-beam fraction
fdir(I)=min(max(exp(3-1.34*Sstar(I)-1.65./max(Sstar(I),1e-6)),0),0.9);
% Solar geometry factors, only where the sun is up (avoids Inf*0 = NaN)
geo_w=zeros(size(Ta)); % wick: tan(theta)/pi
geo_g=zeros(size(Ta)); % globe: 1/(2cos(theta)) - 1
geo_w(I)=tan(theta(I))/pi;
geo_g(I)=1/2./cos(theta(I))-1;

% Transport properties
mu=calcviscosity(18.27e-6,291.15,Ta,120); % dynamic viscosity (kg/m/s)
k=(cp_air+1.25*Rd)*mu; % thermal conductivity of air (W/m/K), Liljegren's thermal_cond

Pr=cp.*mu./k; % Prandtl number
Sc=mu./rho/Df; % Schmidt number
Re=rho.*V*D./mu; % Reynolds number (wick)
h=(k/D)*0.281.*Re.^(1-0.4).*Pr.^(1-0.56); % convective heat transfer coefficient, cylinder (wick)
Reg=rho.*V*Dg./mu; % Reynolds number (globe)
Nu=2.0+0.6*Reg.^0.5.*Pr.^(1/3); % Nusselt number for a sphere, empirically derived
h2=k/Dg.*Nu; % convective heat transfer coefficient derived from Nu (globe)

% Initial guesses
Tw=Ta-(100-RH)/5;
Tg=Ta;

% Calculating Tw and Tg iteratively; each point stops when its own
% change is below tol.  Points with any NaN input are skipped.
valid=~isnan(Ta+RH+S+p+V+theta);
doW=valid; % points where Tw has not converged
doG=valid; % points where Tg has not converged
for iter=1:maxiter
    I=find(doW);
    J=find(doG);
    if isempty(I) && isempty(J)
        break
    end

    % Natural wet bulb (Liljegren eq. 6): net radiation absorbed by the wick
    FA=sigma*eps_w*(0.5*(1+eps_a(I)).*Ta(I).^4-Tw(I).^4)+(1-alpha_w)*S(I).*((1-fdir(I))*(1+lproj)+fdir(I).*(geo_w(I)+lproj)+alpha_sfc);
    ew=611*exp((Lv/Rv)*(1/273.15-1./Tw(I))); % saturation vapor pressure at the wick
    Twnew=Ta(I)-(Lv./cp(I)).*(M_H2O/M_air).*(Pr(I)./Sc(I)).^0.56.*(ew-ea(I))./(p(I)-ew)+FA./h(I);
    Twnew=0.1*Twnew+0.9*Tw(I);
    doW(I)=abs(Tw(I)-Twnew)>=tol;
    Tw(I)=Twnew;

    % Globe temperature (Liljegren eq. 8)
    Tgnew=(0.5*(1+eps_a(J)).*Ta(J).^4-h2(J)/eps_g/sigma.*(Tg(J)-Ta(J))+S(J)/2/eps_g/sigma*(1-alpha_g).*(1+geo_g(J).*fdir(J)+alpha_sfc)).^0.25;
    Tgnew=0.1*Tgnew+0.9*Tg(J);
    doG(J)=abs(Tg(J)-Tgnew)>=tol;
    Tg(J)=Tgnew;
end
if any(doW | doG)
    warning('calcwbgt:noconverge','%d point(s) did not converge after %d iterations; set to NaN.', ...
        nnz(doW | doG),maxiter);
end
Tw(doW | ~valid)=NaN;
Tg(doG | ~valid)=NaN;

% Putting everything together
WBGT = 0.7*Tw + 0.2*Tg + 0.1*Ta - 273.15;
WBGT = reshape(WBGT,sz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q=spechumid(T,RH,p)
% Calculating specific humidity (q) from T (K), RH (percent), p (Pa)
Lv=2.453e6; % Latent heat of vaporization (J/kg)
Rv=461.5; % ideal gas constant for water vapor (J/kg/K)
Rd=287.0; % ideal gas constant for dry air (J/kg/K)
es=611*exp((Lv/Rv)*(1/273.15-1./T));
ws=es*Rd/Rv./(p-es);
w=RH/100.*ws;
q=w./(1+w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mu=calcviscosity(mu0,T0,T,C)
% Dynamic viscosity of air (kg/m/s) from Sutherland's formula;
% C is Sutherlands constant for the gaseous material
mu=mu0*(T0+C)./(T+C).*(T/T0).^1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

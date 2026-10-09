function [zenit,fluxfac]=sza(tau,zday,lat,varargin)

% function [zenit,fluxfac]=sza(tau,zday,lat,[eps,ecc,sol,per],[daypyr,rotper])
%
% This function calculates the cosine of the solar zenith angle based on
% the time of day, the Julian day, and latitude.  You can also specify
% orbital forcings, i.e., obliquity of the Earth's axis of rotation,
% eccentricity of the Earth's orbit, day of the (northern) winter
% solstice, and day of perihelion.  If you don't specify these, the script
% will use present-day values of 23.4393 degrees, 0.0167086, day 355, and
% day 3, respectively.  You can go even farther and specify the number of
% days per year (default is 365) and number of hours per day (default is
% 24).
%
% Ben Kravitz (ben.kravitz.work@gmail.com) 6 February 2012
% Bug fixes 8 October 2026 (with Claude):
%   (1) zenit is now the true cosine of the zenith angle.  It was
%       previously divided by (r/a)^2, i.e., it was cos(zenith) times the
%       Sun-Earth distance flux factor (up to ~3.6% off; could exceed 1).
%       That factor is now returned separately as fluxfac; the old output
%       is zenit.*fluxfac (except where clipped at night).
%   (2) declination now uses sin(delta)=sin(eps)*sin(L), with L the true
%       solar longitude from Kepler's equation (Berger 1978), instead of
%       delta=eps*sin(uniform angle).  Changes delta by up to ~1 degree
%       for the default orbit, more for large eps or ecc.
%   (3) the Sun-Earth distance now uses the eccentric anomaly (exact)
%       instead of a first-order approximation.
%   (4) works for array inputs (^ and * replaced by elementwise operators;
%       tau, zday, and lat may be any compatible sizes).
%   (5) default orbit updated to present-day values (see Defaults):
%       obliquity 23.5 -> 23.4393 deg, eccentricity 0.0178 -> 0.0167086,
%       perihelion day 4 -> 3.  Changes zenit by < 0.002 (declination
%       by < 0.1 deg) and fluxfac by < 0.25%.
%   (6) night values (sun below the horizon) are now exactly 0 instead of
%       1e-8.  Callers that divide by zenit must exclude night points
%       (calcwbgt.m does).
%
% Method
% Mean anomaly M increases uniformly from the day of perihelion; Kepler's
% equation (M = E - ecc*sin(E)) is solved by Newton iteration for the
% eccentric anomaly E, giving the true anomaly nu and r/a = 1 - ecc*cos(E).
% The longitude of perihelion is set so that the true solar longitude
% L = nu + varpi equals 270 degrees on the day of the winter solstice.
% Then sin(delta) = sin(eps)*sin(L), the hour angle is
% H = 2*pi*tau/rotper + pi (local solar time; no equation of time), and
% cos(zenith) = sin(lat)*sin(delta) + cos(lat)*cos(delta)*cos(H).
% The orbital position is evaluated at zday (tau does not advance it).
%
% Defaults (present-day Earth)
% The code takes the days of the winter solstice and perihelion, not the
% longitude of perihelion; the latter follows from them.  Here varpi is
% the Sun's (geocentric) true longitude at perihelion, measured from the
% vernal equinox (= heliocentric longitude of Earth's perihelion + 180).
%   obliquity    23.4393 deg (J2000; 84381.406 arcsec, Capitaine et al.
%                2003, IAU 2006; decreases ~0.013 deg per century)
%   eccentricity 0.0167086 (J2000; Simon et al. 1994)
%   sol = 355    (Dec 21 in a 365-day year; the December solstice falls
%                on Dec 21-22 UTC in the 2020s)
%   per = 3      (Jan 3; perihelion falls on Jan 2-5 in the 2020s)
%   With these, varpi = 283.26 deg and the vernal equinox is day 78.9
%   (Mar 19.9), vs. varpi = 282.94 deg at J2000 (102.94 deg heliocentric,
%   Simon et al. 1994), increasing ~1.72 deg per century (283.4 in 2026).
%   The old defaults (per = 4) gave varpi = 284.3 deg.
%   daypyr = 365 (calendar days, so days 1-365 span one orbit; the
%                tropical year is 365.2422 d)
%   rotper = 24  (mean solar day; tau is local solar time)
%
% References
% Berger, A. L. (1978), Long-term variations of daily insolation and
%   quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
% Capitaine, N., P. T. Wallace, and J. Chapront (2003), Expressions for IAU
%   2000 precession quantities, Astron. Astrophys., 412, 567-586.
% Simon, J. L., P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou,
%   and J. Laskar (1994), Numerical expressions for precession formulae
%   and mean elements for the Moon and the planets, Astron. Astrophys.,
%   282, 663-683.
% Solstice and perihelion dates: U.S. Naval Observatory, "Earth's
%   Seasons and Apsides" (aa.usno.navy.mil/data/EarthSeasons).
% Same geometry as insolation.m, except that insolation.m measures time
% from the vernal equinox and takes the longitude of perihelion as input.
%
% Inputs (tau, zday, lat: scalars or arrays of compatible sizes)
% tau  = time in hours (local solar time; 0 is midnight, 12 is noon, etc.)
% zday = Julian day (day of year, 1 through daypyr)
% lat  = latitude (radians)
% Optional inputs (scalars); if you specify one orbital parameter, you
% have to specify all four:
% eps    = obliquity (radians; default 23.4393*pi/180)
% ecc    = eccentricity (default 0.0167086)
% sol    = Julian day of (northern) winter solstice (default 355)
% per    = Julian day of perihelion (default 3)
% daypyr = number of days per year (default 365)
% rotper = number of hours per day (default 24)
%
% Outputs (size of the implicit expansion of tau, zday, and lat)
% zenit   = cosine of zenith angle.  Night values (cos <= 0) are set to
%           exactly 0, so callers must not divide by zenit at night.
%           NaN inputs give NaN.
% fluxfac = (a/r)^2, the Sun-Earth distance factor (unitless).  Top-of-
%           atmosphere insolation is S0*max(cos(zenith),0)*fluxfac.

nopt=numel(varargin);
if nopt~=0 && nopt~=4 && nopt~=6
    error('sza:nargin',['Specify either no orbital parameters, all four ' ...
        '(eps,ecc,sol,per), or all six (eps,ecc,sol,per,daypyr,rotper).']);
end

% defaults (present-day Earth; see Defaults in the header)
obl=23.4393*pi/180; % obliquity (radians)
ecc=0.0167086;   % eccentricity
sol=355;         % Julian day of winter solstice
per=3;           % Julian day of perihelion
daypyr=365;      % days per year
rotper=24;       % hours per day
if nopt>=4
    obl=varargin{1};
    ecc=varargin{2};
    sol=varargin{3};
    per=varargin{4};
end
if nopt==6
    daypyr=varargin{5};
    rotper=varargin{6};
end
if ~(isscalar(obl) && isscalar(ecc) && isscalar(sol) && isscalar(per) && ...
        isscalar(daypyr) && isscalar(rotper))
    error('sza:scalar','Orbital parameters must be scalars.');
end
if ecc<0 || ecc>=1
    error('sza:ecc','Eccentricity must satisfy 0 <= ecc < 1.');
end

% Orbital position: true anomaly and Sun-Earth distance on day zday
zday=mod(zday,daypyr);
[nu,rdist]=kepler(2*pi*(zday-per)/daypyr,ecc);

% Longitude of perihelion, chosen so that L = 3*pi/2 at the winter solstice
nusol=kepler(2*pi*(sol-per)/daypyr,ecc);
varpi=3*pi/2-nusol;

% Solar declination angle
L=nu+varpi; % true solar longitude (0 at the vernal equinox)
delta=asin(sin(obl)*sin(L));

% Hour angle (H = 2*pi, i.e., cos(H) = 1, at local noon)
tofday=mod(tau,rotper);
H=tofday*2*pi/rotper+pi;

% Calculating zenit = cos(zenith angle); night set to 0
zenit=sin(lat).*sin(delta)+cos(lat).*cos(delta).*cos(H);
zenit(zenit<0)=0; % (unlike max(), keeps NaN inputs as NaN)
fluxfac=1./rdist.^2; % (a/r)^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nu,rdist]=kepler(M,ecc)
% Solves Kepler's equation M = E - ecc*sin(E) by Newton iteration and
% returns the true anomaly nu (radians) and r/a = 1 - ecc*cos(E).
E=M;
for k=1:50
    dE=(E-ecc*sin(E)-M)./(1-ecc*cos(E));
    E=E-dE;
    if all(abs(dE(:))<1e-12)
        break
    end
end
nu=2*atan2(sqrt(1+ecc)*sin(E/2),sqrt(1-ecc)*cos(E/2));
rdist=1-ecc*cos(E);

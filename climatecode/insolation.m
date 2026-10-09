function solarFout=insolation(ecc,eps,oma,lat,day,opts)
% solarFout=insolation(ecc,eps,oma,lat,day)
% solarFout=insolation(ecc,eps,oma,lat,day,Calendar="realworld")
%
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)  12 October 2011
% vectorized by Ben on 12 December 2011
% vectorization fixed by Ben on 16 February 2012
% day 0 changed from the winter solstice to the vernal equinox 8 October 2026
% Calendar option added 8 October 2026 (with Claude); replaces the separate
%   daily_insolation.m (after I. Eisenman), which used the same Berger (1978)
%   equations with the real-world calendar, and the HadCM3-convention
%   insolation.m (Jan 1 = winter solstice, 360-day year)
%
% This function calculates daily mean insolation (W/m2) received at
% the top of the atmosphere at a given latitude and day for specified
% Milankovitch parameters.
% ecc = eccentricity (range of values for Earth: 0.005-0.058)
% eps = obliquity (in radians) (22.1 to 24.5 degrees)
% oma = true solar longitude of perihelion measured from the moving
%       vernal equinox (in radians, 0 to 2pi; present day ~283 degrees)
% lat = latitude (in radians)
% day = day of the year (see Calendar below)
%
% Calendar (optional name-value argument) sets how day is counted:
%   "orbital"   (default) day 0 is the vernal equinox (true solar
%               longitude l=0) and the year is 360 days long.
%   "realworld" day is the calendar day of the year (1 = January 1), the
%               year is 365.2422 days long, and the vernal equinox is fixed
%               at day 80 (21 March).  This is the convention of Eisenman's
%               daily_insolation.m.  Fixing the equinox date is an
%               approximation: for orbits very different from today's the
%               actual equinox date drifts.
%   "hadcm3"    HadCM3 / Unified Model convention: 360-day year with day 1
%               (January 1) at the northern winter solstice (true solar
%               longitude l=270 degrees).
% The solar constant is S0=1365 W/m2; modify the code directly to change it.
%
% Output:  solarFout = daily mean insolation (W/m2), of size
%   [numel(ecc) numel(eps) numel(oma) numel(lat) numel(day)]

%% References
% Berger, Andre L. (1978), Long-term variations of daily
% insolation and quaternary climatic changes, J. Atmos. Sci.,
% 35, 2362-2367,
% doi:10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2.
%
% Peixoto, Jose Pinto and Abraham H. Oort (1992),
% Physics of Climate, 520 pages.

arguments
    ecc {mustBeNumeric,mustBeReal}
    eps {mustBeNumeric,mustBeReal}
    oma {mustBeNumeric,mustBeReal}
    lat {mustBeNumeric,mustBeReal}
    day {mustBeNumeric,mustBeReal}
    opts.Calendar (1,1) string {mustBeMember(opts.Calendar,["orbital","realworld","hadcm3"])} = "orbital"
end

%% parameter list
% lref = true solar longitude at the reference day dayref
switch opts.Calendar
    case "realworld"
        lyear=365.2422; % length of year (in days)
        dayref=80;      % vernal equinox on 21 March
        lref=0;
    case "hadcm3"
        lyear=360.;
        dayref=1;       % winter solstice on January 1
        lref=-pi/2;
    otherwise % "orbital"
        lyear=360.;
        dayref=0;       % vernal equinox on day 0
        lref=0;
end
S0=1365.; % solar constant (in W/m2)

%% vectorization
% output is size [numel(ecc) numel(eps) numel(oma) numel(lat) numel(day)]
[ecc,eps,oma,lat,day]=ndgrid(ecc(:),eps(:),oma(:),lat(:),day(:));

%% convenient quantities to have
ecc2=ecc.^2;
ecc3=ecc.^3;
beta=(1-ecc2).^(0.5);
OMEGA=2*pi/lyear; % angular speed of the earth
                  % also mean longitude angle subtended in one day

%% mean longitude on the reference day (Berger 1978 Section 3)
% lm0 = mean solar longitude where the true longitude is lref (inverse
% equation of center; lref=0 is the vernal equinox, -pi/2 the winter solstice)
lm0term1=(0.5*ecc+0.125*ecc3).*(1+beta).*sin(lref-oma);
lm0term2=(0.25*ecc2).*(0.5+beta).*sin(2*(lref-oma));
lm0term3=(0.125*ecc3).*((1/3)+beta).*sin(3*(lref-oma));
lm0=lref-2*(lm0term1-lm0term2+lm0term3);

%%  This section from Berger 1978 Section 3
% l = solar longitude (angular distance along Earth's orbit)
% l==0 on the vernal equinox
% lm = mean solar longitude (increases uniformly in time from lm0 on the
%      reference day dayref)
% l is then lm plus the equation of center (oma = longitude of perihelion)
dlm=(day-dayref)*OMEGA;
lm=lm0+dlm;
l=lm+(2*ecc-0.25*ecc3).*sin(lm-oma)+1.25*ecc2.*sin(2*(lm-oma))+(13/12)*ecc3.*sin(3*(lm-oma));

%% This section from Peixoto & Oort and Berger 1978 Appendix
% delta = solar declination (lat where sun is directly overhead)
% Ho = hour angle at the horizon (sunrise/sunset)
% rho = Earth-Sun distance in units of the semi-major axis
rho=(1-ecc2)./(1+ecc.*cos(l-oma));
rho2=rho.^2;
delta=asin(sin(eps).*sin(l));
Ho=acos(min(max(-tan(lat).*tan(delta),-1),1)); % pi = polar day, 0 = polar night
solarFout = (S0/pi)*(Ho.*sin(lat).*sin(delta)+cos(lat).*cos(delta).*sin(Ho))./rho2;

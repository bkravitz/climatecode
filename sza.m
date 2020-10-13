function zenit=sza(tau,zday,lat,varargin)

% function zenit=sza(tau,zday,lat,[eps,ecc,sol,per],[daypyr,rotper])
%
% This function calculates the cosine of the solar zenith angle based on
% the time of day, the Julian day, and latitude.  You can also specify
% orbital forcings, i.e., obliquity of the Earth's axis of rotation,
% eccentricity of the Earth's orbit, and day of perihelion.  If you don't
% specify these, the script will use default values of 23.5 degrees,
% 0.0178, and January 1, respectively.  You can go even farther and
% specify the number of days per year (default is 365) and number of
% hours per day (default is 24).
%
% Ben Kravitz (ben.kravitz.work@gmail.com) 6 February 2012
%
% Inputs
% tau  = time in hours (0 is midnight, 12 is noon, etc.)
% zday = Julian day
% lat  = latitude (radians)
% (if number of inputs > 3) eps = obliquity (radians)
% (if number of inputs > 3) ecc = eccentricity
% (if number of inputs > 3) sol = Julian day of winter solstice
% (if number of inputs > 3) per = Julian day of perihelion
% If you want to specify one orbital parameter, you have to specify them all.
% (if number of inputs > 7) daypyr = number of days per year
% (if number of inputs > 7) rotper = number of hours per day
%
% Outputs
% zenit = cosine of zenith angle

if length(varargin)<1;
    eps=23.5*pi/180;
    ecc=0.0178;
    sol=355;
    per=4;
    daypyr=365;
    rotper=24;
elseif length(varargin)<5;
    daypyr=365;
    rotper=24;
    eps=varargin{1};
    ecc=varargin{2};
    sol=varargin{3};
    per=varargin{4};
else
    eps=varargin{1};
    ecc=varargin{2};
    sol=varargin{3};
    per=varargin{4};
    daypyr=varargin{5};
    rotper=varargin{6};
end

% zenit=sin(lat)*sin(delta)+cos(lat)*cos(delta)*cos(H)
% delta=solar declination angle
% H = hour angle

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
rsdist=(1+ecc*cos(2*pi*zdist))^2;

% Calculating zenit
zenit=sin(lat)*sin(delta)/rsdist+cos(lat)*cos(delta)/rsdist*cos(H);
zenit=max(1e-8,zenit);
function output=humidity(T,p,input,flag)

% function output=humidity(T,p,input,flag)
%
% This function converts between relative and specific humidity.
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 13 October 2020
% Bug fixes 8 October 2026 (with Claude):
%   (1) saturation vapor pressure was missing the exp() (it was
%       611.73*(Lv/Rv)*(1/273.16-1/T), not a vapor pressure at all);
%   (2) the specific -> relative humidity branch used an undefined
%       variable q; it is now RH = 100*w/ws;
%   (3) "else if" was missing an end, so the file did not parse.
%   The function could not have run before, so any earlier results did not
%   come from this version.
%
% Method
% Saturation vapor pressure over liquid water from the Clausius-Clapeyron
% equation with constant latent heat,
%   es = 611.73*exp((Lv/Rv)*(1/273.16 - 1/T))   (Pa),
% which differs from Bolton (1980) by -1.3% at 20 C (2306 vs. 2337 Pa),
% -0.7% at 40 C, and +4.5% at -20 C.  Saturation mixing ratio
% ws = (Rd/Rv)*es/(p-es).  Relative humidity is defined as RH = 100*w/ws
% (the mixing-ratio definition; it is smaller than 100*e/es by a factor
% (p-es)/(p-e), e.g., ~2% at 30 C and 50% RH), and q = w/(1+w).
%
% Reference
% Bolton, D. (1980), The computation of equivalent potential temperature,
%   Mon. Wea. Rev., 108, 1046-1053.
%
% Inputs (scalars or arrays of compatible sizes)
% T     = Temperature (K)
% p     = Pressure (Pa)
% input = Humidity (either relative humidity in percent or specific
%         humidity in kg/kg)
% flag  = 1 (Input is relative humidity, output is specific humidity)
%       = 2 (Input is specific humidity, output is relative humidity)
%
% Output
% output = Either specific humidity (kg/kg) or relative humidity (percent),
%          depending on the flag.  Saturation is over liquid water, even
%          below freezing.

if nargin<4 || ~isscalar(flag) || ~any(flag==[1 2])
    error('humidity:flag','flag must be 1 (RH -> q) or 2 (q -> RH).');
end

% Constants
Lv = 2.453e6; % latent heat of vaporization (J/kg)
Rv = 461.5; % ideal gas constant for water vapor (J/kg/K)
Rd = 287.0; % ideal gas constant for dry air (J/kg/K)

es = 611.73*exp((Lv/Rv)*(1/273.16-1./T)); % saturation vapor pressure over water (Pa)
ws = es*(Rd/Rv)./(p-es); % saturation mixing ratio of water (kg/kg)

if flag==1 % relative humidity to specific humidity
    w = (input/100).*ws; % mixing ratio of water from relative humidity
    output = w./(1+w); % specific humidity
else % specific humidity to relative humidity
    w = input./(1-input); % mixing ratio of water from specific humidity
    output = 100*w./ws; % relative humidity
end

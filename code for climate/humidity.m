function output=humidity(T,p,input,flag)

% function output=humidity(T,p,input)
%
% This function converts between relative and specific humidity.
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 13 October 2020
%
% Inputs
% T     = Temperature (K)
% p     = Pressure (Pa)
% input = Humidity (either relative humidity in percent or specific humidity in kg/kg)
% flag  = 1 (Input is relative humidity, output is specific humidity)
%       = 2 (Input is specific humidity, output is relative humidity)
%
% output = Either relative humidity or specific humidity, depending on the flag

% Constants
Lv = 2.453e6; % latent heat of vaporization (J/kg/K)
Rv = 461.5; % ideal gas constant for water vapor (J/kg/K)
Rd = 287.0; % ideal gas constant for dry air (J/kg/K)

es = 611.73*(Lv/Rv)*(1/273.16-1./T); % saturation vapor pressure over water
ws = es*(Rd/Rv)./(p-es); % saturation mixing ratio of water (kg/kg)

if flag==1; % relative humidity to specific humidity
    w = (input/100).*ws; % mixing ratio of water from relative humidity
    output = w./(1+w); % specific humidity
else if flag==2; % specific humidity to relative humidity
    w = input./(1-input); % mixing ratio of water from specific humidity
    output = 100*(1./(1-q))./ws; % relative humidity
end

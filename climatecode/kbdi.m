function Q=kbdi(time,T,rain,R)

% function Q=kbdi(time,T,rain,R)
% Calculates the Keetch-Byram Drought Index
% Coded by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com) on 17 May 2014
% Last modified 10 January 2016
% Bug fixes 8 October 2026 (with Claude): net rainfall now reduces the
%   previous day's index BEFORE the drought factor is computed, so the
%   drought factor uses the post-rain index, as in Keetch and Byram (1968);
%   the drought factor is floored at 0 (previously it went negative for
%   daily maximum temperatures below about 44 F, lowering the index on cold
%   dry days); inputs of different lengths raise an error.  These change Q.
%
% Method: Keetch, J.J. and G.M. Byram (1968), A drought index for forest
% fire control, USDA Forest Service Research Paper SE-38.
%
% Inputs:  time (Julian day plus fractions); can be sub-daily and unsorted,
%               values are grouped into days by floor(time)
%          T (temperature in Fahrenheit); the daily maximum is used
%          rain (in inches); summed over each day
%          R (long-term mean annual precipitation, in inches; scalar)
%          time, T and rain must have the same number of elements.
% Outputs: Q (KBDI; hundredths of an inch of soil moisture drying)
%          Values range between 0 and 800
%          One value per calendar day present in time (sorted, same
%          orientation as unique(floor(time))).  Q starts at 0 (saturated
%          soil) on the first day.  Days with no valid temperature carry
%          the previous value forward (their rain is not used).  NaN rain
%          values are ignored.

if ~isequal(numel(time),numel(T),numel(rain))
    error('kbdi:sizeMismatch','time, T, and rain must have the same number of elements.');
end
[time,J]=sort(time);
T=T(J);
rain=rain(J);
t=unique(floor(time));

Q=zeros(size(t));
Q(1)=0;
rainbucket=0.2; % first 0.2 inch of a rain spell is intercepted and does not reduce Q
for k=2:length(t)
    I=find(floor(time)==t(k));
    TT=max(T(I));
    if isnan(TT)
        Q(k)=Q(k-1);
        continue
    end
    dayrainfall=sum(rain(I),'omitnan');
    if dayrainfall>0
        raintoday=max([dayrainfall-rainbucket 0]); % net rainfall after the 0.2 inch allowance
        rainbucket=max([rainbucket-dayrainfall 0]);
    else
        rainbucket=0.2; % dry day resets the allowance
        raintoday=0;
    end
    Qr=max([Q(k-1)-raintoday*100 0]); % yesterday's index reduced by net rainfall
    dQ=(800-Qr).*(0.968*exp(0.0486*TT)-8.3)./(1+10.88*exp(-0.0441*R))*0.001; % drought factor (Keetch and Byram, 1968)
    Q(k)=Qr+max([dQ 0]); % drought factor cannot be negative
end

function Q=kbdi(time,T,rain,R)

% function Q=kbdi(time,T,rain,R)
% Calculates the Keetch-Byram Drought Index
% Coded by Ben Kravitz (ben.kravitz@pnnl.gov or ben.kravitz.work@gmail.com) on 17 May 2014
% Last modified 10 January 2016
%
% Inputs:  time (Julian day plus fractions)
%          T (temperature in Fahrenheit)
%          rain (in inches)
%          R (long-term mean annual precipitation, in inches)
% Outputs: Q (KBDI; hundredths of an inch of soil moisture drying)
%          Values range between 0 and 800

[time,J]=sort(time);
T=T(J);
rain=rain(J);
t=unique(floor(time));

Q=zeros(size(t));
Q(1)=0;
rainbucket=0.2;
for k=2:length(t);
    I=find(floor(time)==t(k));
    TT=max(T(I));
    if isnan(TT)==1;
        Q(k)=Q(k-1);
        continue;
    end
    if isempty(I)==0;
        dayrainfall=nansum(rain(I));
        if dayrainfall>0;
            raintoday=max([dayrainfall-rainbucket 0]);
            rainbucket=max([rainbucket-dayrainfall 0]);
        else
            rainbucket=0.2;
            raintoday=0;
        end
        dQ=(800-Q(k-1)).*(0.968*exp(0.0486*TT)-8.3)./(1+10.88*exp(-0.0441*R))*0.001; % drought factor
    else
        dQ=0;
    end
    Q(k)=max([Q(k-1)+dQ-raintoday*100 0]);
end

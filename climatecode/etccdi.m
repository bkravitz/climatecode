function indices=etccdi(t,TN,TX,TG,P,TNb,TXb,TGb,Pb,nn,lat)

% function indices=etccdi(t,TN,TX,TG,P,TNb,TXb,TGb,Pb,nn,lat)
%
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 2 February 2020
% Bug fixes 8 August 2026 (with Claude): WSDI/CSDI spell counter reset +
%   end-of-year spell flush; CDD/CWD return 0 (not 1) when no qualifying
%   days; PRCPTOT wet days use >= 1 mm; GSL season-end detection no longer
%   leaves endindex stale (and a season with no qualifying cold spell now
%   runs to the end of the year instead of returning 0); SPII returns 0
%   (not NaN) when a year has no wet days.
% Bug fixes 8 October 2026 (with Claude): R95pTOT/R99pTOT percentiles now use
%   wet days (P >= 1 mm) only; spatial size taken from all but the last dimension.
%   TN10p, TX10p, TN90p, TX90p, WSDI, CSDI now use calendar-day thresholds (5-day
%   window, type 8 percentile) as in ETCCDI/climdex, with the bootstrap applied only
%   to years inside the base period (previously one annual threshold, and every year
%   was bootstrapped).  R95pTOT/R99pTOT use a single wet-day threshold from the full
%   base period (no bootstrap, type 8 percentile).  Years in the percentile indices are
%   now selected with t, like the other indices.
%
% This function takes in a set of daily frequency temperature and precipitation data
% (observed or model output) and produces the 27 ETCCDI metrics defined here:
% http://etccdi.pacificclimate.org/list_27_indices.shtml
%
% The output will be in the form of a matlab structure with similar dimensions to the
% input.  The script is vectorized, so you can input single grid points or entire fields.
% The time dimensions of the outputs will be different, depending on the metrics.
%
% IMPORTANT NOTE:  Please make sure that the time dimension of your inputs is the last
% dimension.  I might add support for arbitrary dimensions later, but I have not done
% that yet.
%
% IMPORTANT NOTE:  This script has year lengths hardcoded to be 365 days.
% I am working on this, but I have not gotten there yet.  For now, if your data has leap
% years, I would recommend dropping February 29 from all applicable years.
% The baseline inputs must have 365-day years.  In the percentile indices (10-15), a
% February 29 in the data is compared against the February 28 threshold.
%
%%%%%%  INPUTS  %%%%%%
% t = time (in Matlab DateVector format), one row per day:  [year month day ...]
% TN = daily minimum temperature (Celsius)
% TX = daily maximum temperature (Celsius)
% TG = daily mean temperature (Celsius)
% P = daily precipitation amount (mm/day)
%     TN, TX, TG, and P are all (space dimensions) x (number of days in t)
% TNb = daily minimum temperature baseline (Celsius) - this must be a complete timeseries, hardcoded to be an integer number of years
% TXb = daily maximum temperature baseline (Celsius) - this must be a complete timeseries, hardcoded to be an integer number of years
% TGb = daily mean temperature baseline (Celsius) - this must be a complete timeseries, hardcoded to be an integer number of years
% Pb = daily precipitation baseline (mm/day) - this must be a complete timeseries, hardcoded to be an integer number of years
%     TNb, TXb, TGb, and Pb are all (space dimensions) x (365*number of base years), starting on Jan 1
%     (TGb is not currently used).  A year of TN (or TX) that is identical to one of the
%     baseline years is treated as inside the base period (see indices 10-15 below), so if
%     the base period overlaps your data, cut the baseline directly out of the data.
% nn = precipitation threshold (mm/day; see #22 below)
% lat = a scalar, vector, or matrix of latitudes (in degrees or radians - does not matter)
%       The size of lat will depend on your input and must be exactly the same size
%       as the non-time component.  For example, if your input is m x n x p with p times,
%       then lat will need to be m x n.  If you are only working with a single grid
%       point, then lat is simply a scalar value - the latitude of that grid point.
%
%%%%%%  OUTPUTS  %%%%%%
% Each field of indices is (space dimensions) x (number of years) for annual indices, or
% (space dimensions) x (12*number of years) for monthly indices (TXX, TNX, TXN, TNN, DTR,
% Rx1day, Rx5day).
% 1:  FD = number of frost days
% 2:  SU = number of summer days
% 3:  ID = number of icing days
% 4:  TR = number of tropical nights
% 5:  GSL = growing season length
% 6:  TXX = monthly maximum value of daily maximum temperature
% 7:  TNX = monthly maximum value of daily minimum temperature
% 8:  TXN = monthly minimum value of daily maximum temperature
% 9:  TNN = monthly minimum value of daily minimum temperature
% 10: TN10p = percentage of days when TN < 10th percentile
% 11: TX10p = percentage of days when TX < 10th percentile
% 12: TN90p = percentage of days when TN > 90th percentile
% 13: TX90p = percentage of days when TX > 90th percentile
% 14: WSDI = warm spell duration index
% 15: CSDI = cold spell duration index
% 16: DTR = daily temperature range
% 17: Rx1day = monthly maximum 1-day precipitation
% 18: Rx5day = monthly maximum consecutive 5-day precipitation
% 19: SPII = simple precipitation intensity index
% 20: R10mm = Annual count of days when precip >= 10 mm
% 21: R20mm = Annual count of days when precip >= 20 mm
% 22: Rnnmm = Annual count of days when precip >= nn mm
% 23: CDD = maximum length of dry spell
% 24: CWD = maximum length of wet spell
% 25: R95pTOT = Annual total precipitation when rain rate > 95th percentile
% 26: R99pTOT = Annual total precipitation when rain rate > 99th percentile
% 27: PRCPTOT = Annual total precipitation in wet days
%
% Percentile thresholds (ETCCDI; Zhang et al. 2005, J. Climate 18:1641; climdex.pcic):
%   Indices 10-15:  for each calendar day, the 10th/90th percentile of the base-period
%     values in a 5-day window centred on that day (5 x number of base years values).
%     Years outside the base period use the thresholds from the full base period.
%     Years inside the base period are bootstrapped:  the year is replaced in the base
%     period by each of the other base years in turn, and the results are averaged.
%   Indices 25-26:  one threshold per grid cell, the 95th/99th percentile of wet days
%     (P >= 1 mm) in the full base period, no bootstrap.
%   All percentiles use the Hyndman & Fan (1996) type 8 estimator, as in climdex.pcic.
% Indices 10-15 are the slow part (the bootstrap is (number of base years - 1) threshold
%   calculations for every year inside the base period).
% For indices 14-15, note that the answers for years inside the base period may not be
%   integers or may be less than 6 due to bootstrapping and subsequent averaging.
%
% GSL, TN10p, TX10p, TN90p, TX90p, WSDI, CSDI, R95pTOT, and R99pTOT do not make sense
% for just one year, so they require more than one year of input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CODE BEGINS HERE      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODIFY AT YOUR OWN RISK  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Doing some calculations I will need for multiple indices
narginchk(11,11);
datasize=size(TN);

% For preallocation, it will be useful to know how many grid cells we will be working
% with.  Time is the last dimension, so the spatial size is everything before it.
if isvector(TN)
    spacesize=1;
else
    spacesize=datasize(1:end-1);
end
if datasize(end)~=size(t,1) && numel(TN)~=size(t,1)
    error('etccdi:size','The last dimension of TN must match the number of rows of t.');
end
if numel(lat)~=prod(spacesize)
    error('etccdi:size','lat must have one value per grid cell (size of the non-time dimensions).');
end

% Which years we cover
allyears=unique(t(:,1));
nyrs=length(allyears);

J=size(TNb);
nbyrs=J(end)/365;
if isvector(TNb)
    nbyrs=numel(TNb)/365;
end
if nbyrs~=round(nbyrs)
    error('etccdi:size','The baseline must be an integer number of 365-day years.');
end

%% Index 1
% FD = number of frost days
% Annual count of days when TN < 0°C

FD=zeros(prod(spacesize),nyrs);
temp=reshape(TN,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    % The following step:
    %   extract all TN values in the selected year
    %   take the floor to get integer values only
    %   reverse the sign so all formerly negative TN values are now at least 1
    %   take the signum so all TN values less than 0°C should end up 1, and all TN values greater than 0°C should end up -1
    %   add 1 and divide by 2, so all TN values less than 0°C should still be 1, and all TN values greater than 0°C should end up 0
    %   then take the floor again to make sure that any TN==0 values do not mess things up
    %   add up the results for all relevant times
    FD(:,k)=sum(floor((sign(-floor(temp(:,I)))+1)/2),2); 
end
FD=reshape(FD,[spacesize nyrs]);
indices.FD=FD;

%% Index 2
% SU = number of summer days
% Annual count of days when TX > 25°C

SU=zeros(prod(spacesize),nyrs);
temp=reshape(TX,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    % The following step:
    %   extract all TX values in the selected year
    %   subtract 25 so all positive values correspond with TX > 25
    %   take the signum so all values are either 1 (TX > 25), -1 (TX<25), or 0 (TX=25)
    %   add 1 and divide by 2, so all values are either 1 (TX > 25), 0 (TX < 25), or 0.5 (TX=25)
    %   take the floor to get rid of TX=25 values
    %   add up the results for all relevant times
    SU(:,k)=sum(floor((sign(temp(:,I)-25)+1)/2),2);
end
SU=reshape(SU,[spacesize nyrs]);
indices.SU=SU;

%% Index 3
% ID = number of icing days
% Annual count of days when TX < 0°C

ID=zeros(prod(spacesize),nyrs);
temp=reshape(TX,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    ID(:,k)=sum(floor((sign(-floor(temp(:,I)))+1)/2),2); 
end
ID=reshape(ID,[spacesize nyrs]);
indices.ID=ID;

%% Index 4
% TR = number of tropical nights
% Annual count of days when TN > 20°C

TR=zeros(prod(spacesize),nyrs);
temp=reshape(TN,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    TR(:,k)=sum(floor((sign(temp(:,I)-20)+1)/2),2);
end
TR=reshape(TR,[spacesize nyrs]);
indices.TR=TR;

%% Index 5
% GSL = growing season length
% Annual (Jan 1 - Dec 31 in NH, Jul 1 - Jun 30 in SH) count between
%   first span of at least 6 days with TG > 5°C and
%   first span of at least 6 days with TG < 5°C in the second half of the year

% every year except final year
GSL=zeros(prod(spacesize),nyrs);
temp=reshape(TG,[prod(spacesize) length(t)]);
for k=1:nyrs-1
    I=find(t(:,1)==allyears(k));
    J=find(t(:,1)==allyears(k+1));
    K=find(t(:,2)<=6);
    L=find(t(:,2)>6);
    for j=1:prod(spacesize)
        if lat(j)>=0
            early_times=intersect(I,K);
            late_times=intersect(I,L);
        else
            early_times=intersect(I,L);
            late_times=intersect(J,K);
        end
        C=(temp(j,early_times)>5);
        D=(temp(j,late_times)>5);
        
        zz=strfind(C,[1 1 1 1 1 1]);
        if isempty(zz)==0
            startindex=zz(1);
        else
            startindex=10000;
        end
        
        zz=strfind(D,[0 0 0 0 0 0]); % first span of 6 consecutive days below 5°C
        if isempty(zz)==1 % no qualifying cold spell: season runs to the end of the year
            endindex=length(D);
        else
            endindex=zz(1)-1; % season ends the day before the cold spell begins
        end
        
        if (startindex>1000) % no growing-season start found
            GSL(j,k)=0;
        else
            GSL(j,k)=endindex+(length(C)-startindex)+1;
        end
    end
end

%% GSL final year
%% NH:  normal
%% SH:  NaN for all values     
I=find(t(:,1)==allyears(nyrs));
K=find(t(:,2)<=6);
L=find(t(:,2)>6);
for j=1:prod(spacesize)
    if lat(j)>=0
        early_times=intersect(I,K);
        late_times=intersect(I,L);
        C=(temp(j,early_times)>5);
        D=(temp(j,late_times)>5);
        
        zz=strfind(C,[1 1 1 1 1 1]);
        if isempty(zz)==0
            startindex=zz(1);
        else
            startindex=10000;
        end
        
        zz=strfind(D,[0 0 0 0 0 0]); % first span of 6 consecutive days below 5°C
        if isempty(zz)==1 % no qualifying cold spell: season runs to the end of the year
            endindex=length(D);
        else
            endindex=zz(1)-1; % season ends the day before the cold spell begins
        end
        
        if (startindex>1000) % no growing-season start found
            GSL(j,nyrs)=0;
        else
            GSL(j,nyrs)=endindex+(length(C)-startindex)+1;
        end
    else
        GSL(j,nyrs)=NaN;
    end
end
GSL=reshape(GSL,[spacesize nyrs]);
indices.GSL=GSL;

%% Index 6
% TXX = monthly maximum value of daily maximum temperature

TXX=zeros(prod(spacesize),nyrs*12);
temp=reshape(TX,[prod(spacesize) length(t)]);
for k=1:nyrs
    for j=1:12
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        TXX(:,(k-1)*12+j)=max(temp(:,intersect(I1,I2)),[],2);
    end
end
TXX=reshape(TXX,[spacesize nyrs*12]);
indices.TXX=TXX;

%% Index 7
% TNX = monthly maximum value of daily minimum temperature

TNX=zeros(prod(spacesize),nyrs*12);
temp=reshape(TN,[prod(spacesize) length(t)]);
for k=1:nyrs
    for j=1:12
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        TNX(:,(k-1)*12+j)=max(temp(:,intersect(I1,I2)),[],2);
    end
end
TNX=reshape(TNX,[spacesize nyrs*12]);
indices.TNX=TNX;

%% Index 8
% TXN = monthly minimum value of daily maximum temperature

TXN=zeros(prod(spacesize),nyrs*12);
temp=reshape(TX,[prod(spacesize) length(t)]);
for k=1:nyrs
    for j=1:12
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        TXN(:,(k-1)*12+j)=min(temp(:,intersect(I1,I2)),[],2);
    end
end
TXN=reshape(TXN,[spacesize nyrs*12]);
indices.TXN=TXN;

%% Index 9
% TNN = monthly minimum value of daily minimum temperature

TNN=zeros(prod(spacesize),nyrs*12);
temp=reshape(TN,[prod(spacesize) length(t)]);
for k=1:nyrs
    for j=1:12
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        TNN(:,(k-1)*12+j)=min(temp(:,intersect(I1,I2)),[],2);
    end
end
TNN=reshape(TNN,[spacesize nyrs*12]);
indices.TNN=TNN;

%% Indices 10-15:  calendar-day percentile thresholds
% Following Zhang et al. (2005, J. Climate 18:1641) and climdex.pcic:
%   - For each calendar day (365-day calendar), the 10th and 90th percentiles are taken
%     from the base-period values in a 5-day window centred on that day, pooled over
%     all base years (5 x nbyrs values).  Windows at the start/end of a base year use
%     the adjacent base year's days (Dec 30 - Jan 2); there is no data before the first
%     or after the last base year, so those few windows have 1-2 fewer values.
%   - Percentiles use the Hyndman & Fan type 8 estimator (as in climdex.pcic).
%   - Years outside the base period are compared against thresholds from the full
%     base period.
%   - Years inside the base period are bootstrapped:  that year is removed from the
%     base period and replaced by each of the other nbyrs-1 base years in turn, the
%     thresholds are recomputed, and the counts are averaged over the nbyrs-1 replicates.
%     This removes the inhomogeneity between in-base and out-of-base years.
%   - A year of data is "in the base period" if it is identical (all grid cells, all 365
%     days) to one of the years of the baseline input; e.g., if TNb was cut out of TN.
%     If no year matches (baseline from a different dataset), there is no bootstrap.
%   - Feb 29 (if present in the data) is compared against the Feb 28 threshold.
%   - Percentages are relative to the number of non-missing days in the year.

% Day of year (1-365) of every data day, by year
cumdays=[0 31 59 90 120 151 181 212 243 273 304 334];
yearidx=cell(nyrs,1);
yeardoy=cell(nyrs,1);
for z=1:nyrs
    I=find(t(:,1)==allyears(z));
    if size(t,2)>=3
        d=cumdays(t(I,2)).'+t(I,3);
        d(t(I,2)==2 & t(I,3)==29)=59; % Feb 29 uses the Feb 28 threshold
    else % no day column: assume consecutive 365-day years starting on Jan 1
        d=(1:numel(I))';
        if numel(I)>365
            error('etccdi:calendar','t needs a day column (DateVector format) for years longer than 365 days.');
        end
    end
    yearidx{z}=I;
    yeardoy{z}=d';
end

%% Indices 10, 12, and 15 (TN)
% TN10p = percentage of days when TN < 10th percentile
% TN90p = percentage of days when TN > 90th percentile
% CSDI = cold spell duration index
% Annual count of days with at least 6 consecutive days when TN < 10th percentile

TN2=reshape(TN,[prod(spacesize) length(t)]);
TNb2=reshape(TNb,[prod(spacesize) nbyrs*365]);
baseyearTN=match_base_years(TN2,TNb2,yearidx,yeardoy,nbyrs);
[TN10p,TN90p,CSDI]=percentile_indices(TN2,TNb2,baseyearTN,yearidx,yeardoy,nbyrs,'lo');
indices.TN10p=reshape(TN10p,[spacesize nyrs]);

%% Indices 11, 13, and 14 (TX)
% TX10p = percentage of days when TX < 10th percentile
% TX90p = percentage of days when TX > 90th percentile
% WSDI = warm spell duration index
% Annual count of days with at least 6 consecutive days when TX > 90th percentile

TX2=reshape(TX,[prod(spacesize) length(t)]);
TXb2=reshape(TXb,[prod(spacesize) nbyrs*365]);
baseyearTX=match_base_years(TX2,TXb2,yearidx,yeardoy,nbyrs);
[TX10p,TX90p,WSDI]=percentile_indices(TX2,TXb2,baseyearTX,yearidx,yeardoy,nbyrs,'hi');
indices.TX10p=reshape(TX10p,[spacesize nyrs]);
indices.TN90p=reshape(TN90p,[spacesize nyrs]);
indices.TX90p=reshape(TX90p,[spacesize nyrs]);
indices.WSDI=reshape(WSDI,[spacesize nyrs]);
indices.CSDI=reshape(CSDI,[spacesize nyrs]);

%% Index 16
% DTR = daily temperature range
% Monthly mean difference between TX and TN

DTR=zeros(prod(spacesize),nyrs*12);
tempTN=reshape(TN,[prod(spacesize) length(t)]);
tempTX=reshape(TX,[prod(spacesize) length(t)]);
for k=1:nyrs
    for j=1:12
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        DTR(:,(k-1)*12+j)=mean(tempTX(:,intersect(I1,I2))-tempTN(:,intersect(I1,I2)),2);
    end
end
DTR=reshape(DTR,[spacesize nyrs*12]);
indices.DTR=DTR;

%% Index 17
% Rx1day = monthly maximum 1 day precipitation

Rx1day=zeros(prod(spacesize),nyrs*12);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    for j=1:12
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        Rx1day(:,(k-1)*12+j)=max(temp(:,intersect(I1,I2)),[],2);
    end
end
Rx1day=reshape(Rx1day,[spacesize nyrs*12]);
indices.Rx1day=Rx1day;

%% Index 18
% Rx5day = monthly maximum consecutive 5-day precipitation

Rx5day=zeros(prod(spacesize),nyrs*12);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    for j=1:12
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        temp2=temp(:,intersect(I1,I2));
        for q=1:prod(spacesize)
            c=0;
            for z=5:size(temp2,2)
                pamt=sum(temp2(q,z-4:z),2);
                c=max([c pamt]);
            end
            Rx5day(q,(k-1)*12+j)=c;
        end
    end
end
Rx5day=reshape(Rx5day,[spacesize nyrs*12]);
indices.Rx5day=Rx5day;

%% Index 19
% SPII = simple precipitation intensity index
% calculated annually, but arbitrary time lengths are fine

SPII=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    temp2=temp(:,I);
    for q=1:prod(spacesize)
        J=temp2(q,:);
        K=find(J>=1);
        if isempty(K)==1 % no wet days: intensity is 0, not 0/0
            SPII(q,k)=0;
        else
            SPII(q,k)=sum(J(K))/length(K);
        end
    end
end
SPII=reshape(SPII,[spacesize nyrs]);
indices.SPII=SPII;

%% Index 20
% R10mm = annual count of days when P >= 10 mm

R10mm=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize)
        J=find(temp(q,I)>=10);
        R10mm(q,k)=length(J);
    end
end
R10mm=reshape(R10mm,[spacesize nyrs]);
indices.R10mm=R10mm;

%% Index 21
% R20mm = annual count of days when P >= 20 mm

R20mm=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize)
        J=find(temp(q,I)>=20);
        R20mm(q,k)=length(J);
    end
end
R20mm=reshape(R20mm,[spacesize nyrs]);
indices.R20mm=R20mm;

%% Index 22
% Rnnmm = annual count of days when P >= nn mm

Rnnmm=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize)
        J=find(temp(q,I)>=nn);
        Rnnmm(q,k)=length(J);
    end
end
Rnnmm=reshape(Rnnmm,[spacesize nyrs]);
indices.Rnnmm=Rnnmm;

%% Index 23
% CDD = maximum length of dry spell
% maximum number of consecutive days with P < 1 mm

CDD=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize)
        J=find(temp(q,I)<1);
        c=1;
        d=1;
        for r=2:length(J)
            if J(r)==J(r-1)+1
                c=c+1;
            else
                d=max([d c]);
                c=1;
            end
        end
        if isempty(J)==1 % no dry days this year
            CDD(q,k)=0;
        else
            CDD(q,k)=max([d c]);
        end
    end
end
CDD=reshape(CDD,[spacesize nyrs]);
indices.CDD=CDD;

%% Index 24
% CWD = maximum length of wet spell
% maximum number of consecutive days with P >= 1 mm

CWD=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize)
        J=find(temp(q,I)>=1);
        c=1;
        d=1;
        for r=2:length(J)
            if J(r)==J(r-1)+1
                c=c+1;
            else
                d=max([d c]);
                c=1;
            end
        end
        if isempty(J)==1 % no wet days this year
            CWD(q,k)=0;
        else
            CWD(q,k)=max([d c]);
        end
    end
end
CWD=reshape(CWD,[spacesize nyrs]);
indices.CWD=CWD;

%% Indices 25 and 26:  wet-day percentile thresholds
% One threshold per grid cell:  the 95th/99th percentile (type 8, as for indices 10-15)
% of all wet days (P >= 1 mm) in the full base period.  No bootstrap (as in ETCCDI).

P2=reshape(P,[prod(spacesize) length(t)]);
Pb2=reshape(Pb,[prod(spacesize) nbyrs*365]);
wet=Pb2';
wet(wet<1)=NaN; % percentile of wet days (P >= 1 mm) only
wet=sort(wet,1); % missing/dry days sort to the end
nwet=sum(~isnan(wet),1);
temp95=quantile_type8(wet,nwet,0.95)';
temp99=quantile_type8(wet,nwet,0.99)';

%% Index 25
% R95pTOT = annual total P when rain rate > 95th percentile
% Only wet days (P >= 1 mm)

R95pTOT=zeros(prod(spacesize),nyrs);
for z=1:nyrs
    outyear=P2(:,yearidx{z});
    outyear(~((outyear>temp95) & (outyear>=1)))=0;
    R95pTOT(:,z)=sum(outyear,2);
end
R95pTOT=reshape(R95pTOT,[spacesize nyrs]);
indices.R95pTOT=R95pTOT;

%% Index 26
% R99pTOT = annual total P when rain rate > 99th percentile
% Only wet days (P >= 1 mm)

R99pTOT=zeros(prod(spacesize),nyrs);
for z=1:nyrs
    outyear=P2(:,yearidx{z});
    outyear(~((outyear>temp99) & (outyear>=1)))=0;
    R99pTOT(:,z)=sum(outyear,2);
end
R99pTOT=reshape(R99pTOT,[spacesize nyrs]);
indices.R99pTOT=R99pTOT;

%% Index 27
% PRCPTOT = Annual total precipitation in wet days (days with P >= 1 mm)

PRCPTOT=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs
    I=find(t(:,1)==allyears(k));
    for j=1:prod(spacesize)
        Q=temp(j,I);
        J=find(Q>=1); % wet days are P >= 1 mm, matching the ETCCDI definition
        PRCPTOT(j,k)=sum(Q(J));
    end
end
PRCPTOT=reshape(PRCPTOT,[spacesize nyrs]);
indices.PRCPTOT=PRCPTOT;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     HELPER FUNCTIONS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function baseyear=match_base_years(X2,Xb2,yearidx,yeardoy,nbyrs)
% For each data year, the baseline year it is identical to (0 if none).
% X2 = data (grid cells x days), Xb2 = baseline (grid cells x nbyrs*365)
nyrs=numel(yearidx);
baseyear=zeros(nyrs,1);
for z=1:nyrs
    d=yeardoy{z};
    keep=[true diff(d)~=0]; % drop Feb 29 (it shares day 59 with Feb 28)
    if ~isequal(d(keep),1:365)
        continue % not a complete year, so it cannot be a baseline year
    end
    I=yearidx{z};
    Xz=X2(:,I(keep));
    for k=1:nbyrs
        if isequaln(Xz,Xb2(:,(k-1)*365+(1:365)))
            baseyear(z)=k;
            break
        end
    end
end
end

function [pctlo,pcthi,spells]=percentile_indices(X2,Xb2,baseyear,yearidx,yeardoy,nbyrs,spellside)
% Percentage of days below the 10th / above the 90th calendar-day percentile, and the
% number of days in spells of at least 6 days below the 10th ('lo') or above the
% 90th ('hi') percentile, for every grid cell and year.  See the notes above index 10.
ngrid=size(X2,1);
nyrs=numel(yearidx);
pctlo=zeros(ngrid,nyrs);
pcthi=zeros(ngrid,nyrs);
spells=zeros(ngrid,nyrs);
% Work on blocks of grid cells so each window array is roughly 80 MB
chunk=max(1,floor(1e7/(5*nbyrs*365)));
for g1=1:chunk:ngrid
    G=g1:min(g1+chunk-1,ngrid);
    W=window_samples(Xb2(G,:),nbyrs);
    Ws=sort(W,1); % each base year's 5 window values, sorted (for the bootstrap)
    nWs=sum(~isnan(Ws),1); % number of non-missing values in each
    Ws(isnan(Ws))=Inf; % missing values last, keeping the columns sorted
    [q10,q90]=window_thresholds(W,1:nbyrs,numel(G)); % full base period
    for z=1:nyrs
        x=X2(G,yearidx{z});
        d=yeardoy{z};
        k=baseyear(z);
        if k==0 || nbyrs<2 % out of base:  full-base thresholds
            [nlo,nhi,nsp]=count_exceedances(x,q10(:,d),q90(:,d),spellside);
        else % in base:  replace year k by each other base year j, then average
            % Each bootstrap sample is (base period without year k) + (year j again),
            % so sort the base period without year k once and merge in year j each time.
            others=setdiff(1:nbyrs,k);
            A=sort(reshape(W(:,others,:),5*(nbyrs-1),[]),1);
            nA=sum(~isnan(A),1);
            A(isnan(A))=Inf;
            nlo=0;
            nhi=0;
            nsp=0;
            for j=others
                [r10,r90]=merged_thresholds(A,reshape(Ws(:,j,:),5,[]),nA+reshape(nWs(1,j,:),1,[]),numel(G));
                [a,b,c]=count_exceedances(x,r10(:,d),r90(:,d),spellside);
                nlo=nlo+a;
                nhi=nhi+b;
                nsp=nsp+c;
            end
            nlo=nlo/(nbyrs-1);
            nhi=nhi/(nbyrs-1);
            nsp=nsp/(nbyrs-1);
        end
        ndays=sum(~isnan(x),2);
        pctlo(G,z)=100*nlo./ndays;
        pcthi(G,z)=100*nhi./ndays;
        spells(G,z)=nsp;
    end
end
end

function W=window_samples(Xb2,nbyrs)
% Base-period values in the 5-day window around each calendar day, for each base year.
% Output is 5 x nbyrs x (grid cells*365), with grid cells varying fastest.
% The baseline is treated as one continuous series, so windows cross year boundaries;
% days before the first or after the last base year are NaN.
ngrid=size(Xb2,1);
Bp=[NaN(ngrid,2) Xb2 NaN(ngrid,2)];
[o,y,d]=ndgrid(-2:2,1:nbyrs,1:365);
col=(y-1)*365+d+o+2; % column in the padded series
W=reshape(Bp(:,col(:)),[ngrid 5 nbyrs 365]);
W=reshape(permute(W,[2 3 1 4]),[5 nbyrs ngrid*365]);
end

function [q10,q90]=window_thresholds(W,yrs,ngrid)
% 10th and 90th percentiles (grid cells x 365) of the window samples from base years yrs
X=sort(reshape(W(:,yrs,:),5*numel(yrs),[]),1); % missing values sort to the end
n=sum(~isnan(X),1);
q10=reshape(quantile_type8(X,n,0.1),[ngrid 365]);
q90=reshape(quantile_type8(X,n,0.9),[ngrid 365]);
end

function [q10,q90]=merged_thresholds(A,b,n,ngrid)
% Same as window_thresholds for the sample [A; b] without sorting it:  A (M x columns)
% and b (5 x columns) are each sorted along columns, with missing values set to Inf,
% and n is the number of non-missing values in each column of [A; b].
% Uses r-th smallest of [A; b] = min over i=0..5 of max(A(r-i),b(i)), A(0)=b(0)=-Inf.
M=size(A,1);
offset=(0:size(A,2)-1)*M;
q=cell(1,2);
p=[0.1 0.9];
for m=1:2
    [lo,hi,g]=type8_positions(n,p(m));
    xr=cell(1,2);
    r=[lo;hi];
    for s=1:2
        xr{s}=Inf(size(n));
        for i=0:5
            ia=r(s,:)-i; % take ia values from A and i values from b
            a=A(min(max(ia,1),M)+offset);
            a(ia<1)=-Inf;
            if i>0
                a=max(a,b(i,:));
            end
            a(ia<0 | ia>M)=Inf; % impossible split
            xr{s}=min(xr{s},a);
        end
    end
    q{m}=xr{1}+g.*(xr{2}-xr{1});
    q{m}(n==0)=NaN;
end
q10=reshape(q{1},[ngrid 365]);
q90=reshape(q{2},[ngrid 365]);
end

function q=quantile_type8(xs,n,p)
% Hyndman & Fan (1996) type 8 quantile (median-unbiased; used by climdex.pcic) of each
% column of xs, which must be sorted ascending with the n non-missing values first.
% Same as R's quantile(x,p,type=8).
[lo,hi,g]=type8_positions(n,p);
offset=(0:size(xs,2)-1)*size(xs,1);
q=xs(lo+offset)+g.*(xs(hi+offset)-xs(lo+offset));
q(n==0)=NaN;
end

function [lo,hi,g]=type8_positions(n,p)
% Type 8 quantile p of n sorted values is x(lo)+g*(x(hi)-x(lo))
fuzz=4*eps;
h=n*p+(p+1)/3;
j=floor(h+fuzz);
g=h-j;
g(abs(g)<fuzz)=0;
lo=min(max(j,1),max(n,1));
hi=min(max(j+1,1),max(n,1));
end

function [nlo,nhi,nsp]=count_exceedances(x,thlo,thhi,spellside)
% Days below thlo, days above thhi, and days in spells of at least 6 consecutive
% exceedances (of thlo if spellside is 'lo', of thhi if 'hi')
below=x<thlo;
above=x>thhi;
nlo=sum(below,2);
nhi=sum(above,2);
if strcmp(spellside,'lo')
    nsp=spell_days(below);
else
    nsp=spell_days(above);
end
end

function n=spell_days(E)
% Number of true values in each row of E that belong to a run of at least 6
n=zeros(size(E,1),1);
if size(E,2)<6
    return
end
start6=movsum(double(E),[0 5],2,'Endpoints','discard')==6; % days s to s+5 all true
inspell=movmax(double([start6 false(size(E,1),5)]),[5 0],2); % day lies in such a run
n=sum(inspell,2);
end

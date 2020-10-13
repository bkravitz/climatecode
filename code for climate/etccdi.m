function indices=etccdi(t,TN,TX,TG,P,TNb,TXb,TGb,Pb,nn,lat)

% function indices=etccdi(t,TN,TX,TG,P,TNb,TXb,TGb,Pb,nn)
%
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 2 February 2020
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
%
%%%%%%  INPUTS  %%%%%%
% t = time (in Matlab DateVector format)
% TN = daily minimum temperature (Celsius)
% TX = daily maximum temperature (Celsius)
% TG = daily mean temperature (Celsius)
% P = daily precipitation amount (mm/day)
% TNb = daily minimum temperature baseline (Celsius) - this must be a complete timeseries, hardcoded to be an integer number of years
% TXb = daily maximum temperature baseline (Celsius) - this must be a complete timeseries, hardcoded to be an integer number of years
% TGb = daily mean temperature baseline (Celsius) - this must be a complete timeseries, hardcoded to be an integer number of years
% Pb = daily precipitation baseline (mm/day) - this must be a complete timeseries, hardcoded to be an integer number of years
% nn = precipitation threshold (mm/day; see #22 below)
% lat = a scalar, vector, or matrix of latitudes (in degrees or radians - does not matter)
%       The size of lat will depend on your input and must be exactly the same size
%       as the non-time component.  For example, if your input is m x n x p with p times,
%       then lat will need to be m x n.  If you are only working with a single grid
%       point, then lat is simply a scalar value - the latitude of that grid point.
%
%%%%%%  OUTPUTS  %%%%%%
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
% 20: R10mm = Annual count of days when precip > 10 mm
% 21: R20mm = Annual count of days when precip > 20 mm
% 22: Rnnmm = Annual count of days when precip > nn mm
% 23: CDD = maximum length of dry spell
% 24: CWD = maximum length of wet spell
% 25: R95pTOT = Annual total precipitation when rain rate > 95th percentile
% 26: R99pTOT = Annual total precipitation when rain rate > 99th percentile
% 27: PRCPTOT = Annual total precipitation in wet days
%
% Note that indices 10-15 and 25-26 are slow (approx. 0.4 seconds per grid cell per year each)
% For indices 14-15, note that the answers may not be integers or may be less than 6
%   due to bootstrapping and subsequent averaging.
%
% GSL, TN10p, TX10p, TN90p, TX90p, WSDI, CSDI, R95pTOT, and R99pTOT do not make sense
% for just one year, so they require more than one year of input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CODE BEGINS HERE      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODIFY AT YOUR OWN RISK  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Doing some calculations I will need for multiple indices
datasize=size(TN);
basesize=size(TNb);

% For preallocation, it will be useful to know how many grid cells we will be working
% with.  This requires figuring out which is the time dimension and not including it.
n=ndims(TN);
if n>1;
    timedim=find(datasize==length(t));
    spacedim=setdiff(1:n,timedim);
    spacesize=datasize(spacedim);
else
    spacesize=1;
end

% Which years we cover
allyears=unique(t(:,1));
nyrs=length(allyears);

J=size(TNb);
nbyrs=J(end)/365;

molens=[31 28 31 30 31 30 31 31 30 31 30 31]; % hardcoding no leap days

%% Index 1
% FD = number of frost days
% Annual count of days when TN < 0°C

FD=zeros(prod(spacesize),nyrs);
temp=reshape(TN,[prod(spacesize) length(t)]);
for k=1:nyrs;
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
for k=1:nyrs;
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
for k=1:nyrs;
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
for k=1:nyrs;
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
for k=1:nyrs-1;
    I=find(t(:,1)==allyears(k));
    J=find(t(:,1)==allyears(k+1));
    K=find(t(:,2)<=6);
    L=find(t(:,2)>6);
    for j=1:prod(spacesize);
        if lat(j)>=0;
            early_times=intersect(I,K);
            late_times=intersect(I,L);
        else
            early_times=intersect(I,L);
            late_times=intersect(J,K);
        end
        C=(temp(j,early_times)>5);
        D=(temp(j,late_times)>5);
        
        zz=strfind(C,[1 1 1 1 1 1]);
        if isempty(zz)==0;
            startindex=zz(1);
        else
            startindex=10000;
        end
        
        zz=strfind(D,[1 0]); % looking for transitions from above to below 5°C
        if isempty(zz)==1; % case with no transitions
            if sum(D)==length(D); % everything is above 5°C
                endindex=length(D);
            elseif sum(D)==0; % everything is below 5°C
                endindex=10000;
            end
        else % there is a transition
            zz=strfind(D,[0 0 0 0 0 0]); % finding spans of 6 days
            if isempty(zz)==1; % did not find anything
                endindex=length(D);
            else
                endindex=zz(1)-1;
            end
        end
        
        if ((startindex>1000) | (endindex>1000));
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
for j=1:prod(spacesize);
    if lat(j)>=0;
        early_times=intersect(I,K);
        late_times=intersect(I,L);
        C=(temp(j,early_times)>5);
        D=(temp(j,late_times)>5);
        
        zz=strfind(C,[1 1 1 1 1 1]);
        if isempty(zz)==0;
            startindex=zz(1);
        else
            startindex=10000;
        end
        
        zz=strfind(D,[1 0]); % looking for transitions from above to below 5°C
        if isempty(zz)==1; % case with no transitions
            if sum(D)==length(D); % everything is above 5°C
                endindex=length(D);
            elseif sum(D)==0; % everything is below 5°C
                endindex=10000;
            end
        else % there is a transition
            zz=strfind(D,[0 0 0 0 0 0]); % finding spans of 6 days
            if isempty(zz)==1; % did not find anything
                endindex=length(D);
            else
                endindex=zz(1)-1;
            end
        end
        
        if ((startindex>1000) | (endindex>1000));
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
for k=1:nyrs;
    for j=1:12;
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
for k=1:nyrs;
    for j=1:12;
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
for k=1:nyrs;
    for j=1:12;
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        TXN(:,(k-1)*12+j)=min(temp(:,intersect(I1,I2)),[],2);
    end
end
TXN=reshape(TXN,[spacesize nyrs*12]);
indices.TXN=TXN;

%% Index 9
% TNN = monthly minimum value of daily maximum temperature

TNN=zeros(prod(spacesize),nyrs*12);
temp=reshape(TN,[prod(spacesize) length(t)]);
for k=1:nyrs;
    for j=1:12;
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        TNN(:,(k-1)*12+j)=min(temp(:,intersect(I1,I2)),[],2);
    end
end
TNN=reshape(TNN,[spacesize nyrs*12]);
indices.TNN=TNN;

%% Index 10
% TN10p = percentage of days when TN < 10th percentile

TN2=reshape(TN,[prod(spacesize) length(t)]);
TNb2=reshape(TNb,[prod(spacesize) nbyrs*365]);
TN10p=zeros(prod(spacesize),nyrs);
for z=1:nyrs;
    outyear=TN2(:,(z-1)*365+1:(z-1)*365+365);
    exceedance10=zeros(prod(spacesize),nbyrs);
    for k=1:nbyrs;
        outtimes=(k-1)*365+1:(k-1)*365+365;
        keeptimes=setdiff(1:365*nbyrs,outtimes);
        keepyears=TNb2(:,keeptimes);
        exceedtemp=zeros(prod(spacesize),nbyrs-1);
        for j=1:(nbyrs-1);
            keepyears2=cat(2,keepyears,keepyears(:,(j-1)*365+1:(j-1)*365+365));
            temp10=quantile(keepyears2,0.1,2);
            for q=1:prod(spacesize);
                I=find(outyear(q,:)<temp10(q));
                exceedtemp(q,j)=length(I);
            end
        end
        exceedance10(:,k)=mean(exceedtemp,2);
    end
    TN10p(:,z)=mean(exceedance10,2);
end
TN10p=reshape(100*TN10p/365,[spacesize nyrs]);
indices.TN10p=TN10p;

%% Index 11
% TX10p = percentage of days when TX < 10th percentile

TX2=reshape(TX,[prod(spacesize) length(t)]);
TXb2=reshape(TXb,[prod(spacesize) nbyrs*365]);
TX10p=zeros(prod(spacesize),nyrs);
for z=1:nyrs;
    outyear=TX2(:,(z-1)*365+1:(z-1)*365+365);
    exceedance10=zeros(prod(spacesize),nbyrs);
    for k=1:nbyrs;
        outtimes=(k-1)*365+1:(k-1)*365+365;
        keeptimes=setdiff(1:365*nbyrs,outtimes);
        keepyears=TXb2(:,keeptimes);
        exceedtemp=zeros(prod(spacesize),nbyrs-1);
        for j=1:(nbyrs-1);
            keepyears2=cat(2,keepyears,keepyears(:,(j-1)*365+1:(j-1)*365+365));
            temp10=quantile(keepyears2,0.1,2);
            for q=1:prod(spacesize);
                I=find(outyear(q,:)<temp10(q));
                exceedtemp(q,j)=length(I);
            end
        end
        exceedance10(:,k)=mean(exceedtemp,2);
    end
    TX10p(:,z)=mean(exceedance10,2);
end
TX10p=reshape(100*TX10p/365,[spacesize nyrs]);
indices.TX10p=TX10p;

%% Index 12
% TN90p = percentage of days when TN > 90th percentile

TN2=reshape(TN,[prod(spacesize) length(t)]);
TNb2=reshape(TNb,[prod(spacesize) nbyrs*365]);
TN90p=zeros(prod(spacesize),nyrs);
for z=1:nyrs;
    outyear=TN2(:,(z-1)*365+1:(z-1)*365+365);
    exceedance90=zeros(prod(spacesize),nbyrs);
    for k=1:nbyrs;
        outtimes=(k-1)*365+1:(k-1)*365+365;
        keeptimes=setdiff(1:365*nbyrs,outtimes);
        keepyears=TNb2(:,keeptimes);
        exceedtemp=zeros(prod(spacesize),nbyrs-1);
        for j=1:(nbyrs-1);
            keepyears2=cat(2,keepyears,keepyears(:,(j-1)*365+1:(j-1)*365+365));
            temp90=quantile(keepyears2,0.9,2);
            for q=1:prod(spacesize);
                I=find(outyear(q,:)>temp90(q));
                exceedtemp(q,j)=length(I);
            end
        end
        exceedance90(:,k)=mean(exceedtemp,2);
    end
    TN90p(:,z)=mean(exceedance90,2);
end
TN90p=reshape(100*TN90p/365,[spacesize nyrs]);
indices.TN90p=TN90p;

%% Index 13
% TX90p = percentage of days when TX > 90th percentile

TX2=reshape(TX,[prod(spacesize) length(t)]);
TXb2=reshape(TXb,[prod(spacesize) nbyrs*365]);
TX90p=zeros(prod(spacesize),nyrs);
for z=1:nyrs;
    outyear=TX2(:,(z-1)*365+1:(z-1)*365+365);
    exceedance90=zeros(prod(spacesize),nbyrs);
    for k=1:nbyrs;
        outtimes=(k-1)*365+1:(k-1)*365+365;
        keeptimes=setdiff(1:365*nbyrs,outtimes);
        keepyears=TXb2(:,keeptimes);
        exceedtemp=zeros(prod(spacesize),nbyrs-1);
        for j=1:(nbyrs-1);
            keepyears2=cat(2,keepyears,keepyears(:,(j-1)*365+1:(j-1)*365+365));
            temp90=quantile(keepyears2,0.9,2);
            for q=1:prod(spacesize);
                I=find(outyear(q,:)>temp90(q));
                exceedtemp(q,j)=length(I);
            end
        end
        exceedance90(:,k)=mean(exceedtemp,2);
    end
    TX90p(:,z)=mean(exceedance90,2);
end
TX90p=reshape(100*TX90p/365,[spacesize nyrs]);
indices.TX90p=TX90p;

%% Index 14
% WSDI = warm spell duration index
% Annual count of days with at least 6 consecutive days when TX > 90th percentile

TX2=reshape(TX,[prod(spacesize) length(t)]);
TXb2=reshape(TXb,[prod(spacesize) nbyrs*365]);
WSDI=zeros(prod(spacesize),nyrs);
for z=1:nyrs;
    outyear=TX2(:,(z-1)*365+1:(z-1)*365+365);
    exceedance90=zeros(prod(spacesize),nbyrs);
    for k=1:nbyrs;
        outtimes=(k-1)*365+1:(k-1)*365+365;
        keeptimes=setdiff(1:365*nbyrs,outtimes);
        keepyears=TXb2(:,keeptimes);
        exceedtemp=zeros(prod(spacesize),nbyrs-1);
        for j=1:(nbyrs-1);
            keepyears2=cat(2,keepyears,keepyears(:,(j-1)*365+1:(j-1)*365+365));
            temp90=quantile(keepyears2,0.9,2);
            for q=1:prod(spacesize);
                I=find(outyear(q,:)>temp90(q));
                c=1;
                for r=2:length(I);
                    if I(r)==I(r-1)+1;
                        c=c+1;
                    else
                        if c>5;
                            exceedtemp(q,j)=exceedtemp(q,j)+c;
                            c=1;
                        end
                    end
                end
            end
        end
        exceedance90(:,k)=mean(exceedtemp,2);
    end
    WSDI(:,z)=mean(exceedance90,2);
end
WSDI=reshape(WSDI,[spacesize nyrs]);
indices.WSDI=WSDI;

%% Index 15
% CSDI = cold spell duration index
% Annual count of days with at least 6 consecutive days when TN < 10th percentile

TN2=reshape(TN,[prod(spacesize) length(t)]);
TNb2=reshape(TNb,[prod(spacesize) nbyrs*365]);
CSDI=zeros(prod(spacesize),nyrs);
for z=1:nyrs;
    outyear=TN2(:,(z-1)*365+1:(z-1)*365+365);
    exceedance10=zeros(prod(spacesize),nbyrs);
    for k=1:nbyrs;
        outtimes=(k-1)*365+1:(k-1)*365+365;
        keeptimes=setdiff(1:365*nbyrs,outtimes);
        keepyears=TNb2(:,keeptimes);
        exceedtemp=zeros(prod(spacesize),nbyrs-1);
        for j=1:(nbyrs-1);
            keepyears2=cat(2,keepyears,keepyears(:,(j-1)*365+1:(j-1)*365+365));
            temp10=quantile(keepyears2,0.1,2);
            for q=1:prod(spacesize);
                I=find(outyear(q,:)<temp10(q));
                c=1;
                for r=2:length(I);
                    if I(r)==I(r-1)+1;
                        c=c+1;
                    else
                        if c>5;
                            exceedtemp(q,j)=exceedtemp(q,j)+c;
                            c=1;
                        end
                    end
                end
            end
        end
        exceedance10(:,k)=mean(exceedtemp,2);
    end
    CSDI(:,z)=mean(exceedance10,2);
end
CSDI=reshape(CSDI,[spacesize nyrs]);
indices.CSDI=CSDI;

%% Index 16
% DTR = daily temperature range
% Monthly mean difference between TX and TN

DTR=zeros(prod(spacesize),nyrs*12);
tempTN=reshape(TN,[prod(spacesize) length(t)]);
tempTX=reshape(TX,[prod(spacesize) length(t)]);
for k=1:nyrs;
    for j=1:12;
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
for k=1:nyrs;
    for j=1:12;
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
for k=1:nyrs;
    for j=1:12;
        I1=find(t(:,1)==allyears(k));
        I2=find(t(:,2)==j);
        temp2=temp(:,intersect(I1,I2));
        for q=1:prod(spacesize);
            c=0;
            for z=5:molens(j);
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
for k=1:nyrs;
    I=find(t(:,1)==allyears(k));
    temp2=temp(:,I);
    for q=1:prod(spacesize);
        J=temp2(q,:);
        K=find(J>=1);
        SPII(q,k)=sum(J(K))/length(K);
    end
end
SPII=reshape(SPII,[spacesize nyrs]);
indices.SPII=SPII;

%% Index 20
% R10mm = annual count of days when P >= 10 mm

R10mm=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs;
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize);
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
for k=1:nyrs;
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize);
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
for k=1:nyrs;
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize);
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
for k=1:nyrs;
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize);
        J=find(temp(q,I)<1);
        c=1;
        d=1;
        for r=2:length(J);
            if J(r)==J(r-1)+1;
                c=c+1;
            else
                d=max([d c]);
                c=1;
            end
        end
        CDD(q,k)=max([d c]);
    end
end
CDD=reshape(CDD,[spacesize nyrs]);
indices.CDD=CDD;

%% Index 24
% CWD = maximum length of wet spell
% maximum number of consecutive days with P >= 1 mm

CWD=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs;
    I=find(t(:,1)==allyears(k));
    for q=1:prod(spacesize);
        J=find(temp(q,I)>=1);
        c=1;
        d=1;
        for r=2:length(J);
            if J(r)==J(r-1)+1;
                c=c+1;
            else
                d=max([d c]);
                c=1;
            end
        end
        CWD(q,k)=max([d c]);
    end
end
CWD=reshape(CWD,[spacesize nyrs]);
indices.CWD=CWD;

%% Index 25
% R95pTOT = annual total P when rain rate > 95th percentile
% Only wet days (P >= 1 mm)

P2=reshape(P,[prod(spacesize) length(t)]);
Pb2=reshape(Pb,[prod(spacesize) nbyrs*365]);
R95pTOT=zeros(prod(spacesize),nyrs);
for z=1:nyrs;
    outyear=P2(:,(z-1)*365+1:(z-1)*365+365);
    exceedance95=zeros(prod(spacesize),nbyrs);
    for k=1:nbyrs;
        outtimes=(k-1)*365+1:(k-1)*365+365;
        keeptimes=setdiff(1:365*nbyrs,outtimes);
        keepyears=Pb2(:,keeptimes);
        exceedtemp=zeros(prod(spacesize),nbyrs-1);
        for j=1:(nbyrs-1);
            keepyears2=cat(2,keepyears,keepyears(:,(j-1)*365+1:(j-1)*365+365));
            temp95=quantile(keepyears2,0.95,2);
            for q=1:prod(spacesize);
                I=find((outyear(q,:)>temp95(q)) & (outyear(q,:)>=1));
                exceedtemp(q,j)=sum(outyear(q,I),2);
            end
        end
        exceedance95(:,k)=mean(exceedtemp,2);
    end
    R95pTOT(:,z)=mean(exceedance95,2);
end
R95pTOT=reshape(R95pTOT,[spacesize nyrs]);
indices.R95pTOT=R95pTOT;

%% Index 26
% R99pTOT = annual total P when rain rate > 99th percentile
% Only wet days (P >= 1 mm)

P2=reshape(P,[prod(spacesize) length(t)]);
Pb2=reshape(Pb,[prod(spacesize) nbyrs*365]);
R99pTOT=zeros(prod(spacesize),nyrs);
for z=1:nyrs;
    outyear=P2(:,(z-1)*365+1:(z-1)*365+365);
    exceedance99=zeros(prod(spacesize),nbyrs);
    for k=1:nbyrs;
        outtimes=(k-1)*365+1:(k-1)*365+365;
        keeptimes=setdiff(1:365*nbyrs,outtimes);
        keepyears=Pb2(:,keeptimes);
        exceedtemp=zeros(prod(spacesize),nbyrs-1);
        for j=1:(nbyrs-1);
            keepyears2=cat(2,keepyears,keepyears(:,(j-1)*365+1:(j-1)*365+365));
            temp99=quantile(keepyears2,0.99,2);
            for q=1:prod(spacesize);
                I=find((outyear(q,:)>temp99(q)) & (outyear(q,:)>=1));
                exceedtemp(q,j)=sum(outyear(q,I),2);
            end
        end
        exceedance99(:,k)=mean(exceedtemp,2);
    end
    R99pTOT(:,z)=mean(exceedance99,2);
end
R99pTOT=reshape(R99pTOT,[spacesize nyrs]);
indices.R99pTOT=R99pTOT;

%% Index 27
% PRCPTOT = Annual total precipitation in wet days (days with P >= 1 mm)

PRCPTOT=zeros(prod(spacesize),nyrs);
temp=reshape(P,[prod(spacesize) length(t)]);
for k=1:nyrs;
    I=find(t(:,1)==allyears(k));
    for j=1:prod(spacesize);
        Q=temp(j,I);
        J=find(Q>1);
        PRCPTOT(j,k)=sum(Q(J));
    end
end
PRCPTOT=reshape(PRCPTOT,[spacesize nyrs]);
indices.PRCPTOT=PRCPTOT;

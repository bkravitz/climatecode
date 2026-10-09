function [r,n,rmod,rmed,rbar,stdev,reff]=makelogndistro(qt1,qt2,flag,rmin,rmax,nr)
% [r,n,rmod,rmed,rbar,stdev,reff]=makelogndistro(qt1,qt2,flag,rmin,rmax,nr)
%
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)  12 October 2011
% Bug fixes 8 October 2026 (with Claude): for flag=3 the mean radius is now
%   solved exactly with fzero (previously it was taken from the radius grid
%   where the residual changed sign, which was only as accurate as the grid
%   and could return several values or none); n is computed directly from
%   the lognormal formula (no Statistics Toolbox needed; same values);
%   inputs are validated.  Only flag=3 results change.  Flags 1, 2 and 4
%   now raise an error unless qt2>qt1 (previously qt2<qt1 gave a complex
%   sigma and qt2==qt1 gave sigma=0 and a NaN distribution).
%
% generates a lognormal distribution with inputs:
% - qt1 and qt2 are any two quantities from the list:
%     mode radius, median radius, mean radius, standard deviation
%     (or the combination of geometric mean and standard deviation)
%     of the lognormal distribution; radii in any consistent units
%     (e.g. microns), and all radii outputs are in the same units
% - flag specifies what qt1 and qt2 are:
%     1 - qt1=mode radius,    qt2=median radius
%     2 - qt1=mode radius,    qt2=mean radius
%     3 - qt1=mode radius,    qt2=standard deviation
%     4 - qt1=median radius,  qt2=mean radius
%     5 - qt1=median radius,  qt2=standard deviation
%     6 - qt1=mean radius,    qt2=standard deviation
%     7 - qt1=geometric mean, qt2=geometric standard deviation
%     flag=7 are commonly used quantities in reporting size distributions
%     (flags 1, 2 and 4 require qt2>qt1 strictly, since mode < median <
%     mean for any lognormal with sigma>0; otherwise an error is raised)
% - rmin and rmax (limits of calculation of the distribution)
%     I recommend 1e-4 and 1e4, respectively
% - nr (number of radii in the interval [rmin,rmax])
%     I recommend 1e6
%
% outputs:
% - r (list of radii, 1 x nr, logarithmically spaced)
% - n (size distribution dN/dr at each r, 1 x nr; the lognormal pdf in r,
%     which integrates to 1 over 0<r<inf, in units of 1/radius)
% - rmod (mode radius)
% - rmed (median radius, also geometric mean)
% - rbar (mean radius)
% - stdev (standard deviation of the radius, not the geometric one)
% - reff (aerosol effective radius, sum(r^3*n*dr)/sum(r^2*n*dr), computed
%     numerically over [rmin,rmax]; Hansen and Travis, 1974)
%
% Lognormal relations used: with mu=log(rmed) and sigma=log(geometric
% standard deviation), rmod=exp(mu-sigma^2), rbar=exp(mu+sigma^2/2),
% stdev^2=rbar^2*(exp(sigma^2)-1).

arguments
    qt1 (1,1) {mustBeNumeric,mustBePositive}
    qt2 (1,1) {mustBeNumeric,mustBePositive}
    flag (1,1) {mustBeMember(flag,1:7)}
    rmin (1,1) {mustBeNumeric,mustBePositive}
    rmax (1,1) {mustBeNumeric,mustBeGreaterThan(rmax,rmin)}
    nr (1,1) {mustBeInteger,mustBeGreaterThan(nr,1)}
end

% flags 1, 2 and 4 take sigma^2 proportional to log(qt2/qt1), which must be > 0
if ismember(flag,[1 2 4]) && ~(qt2>qt1)
    names={'mode radius','median radius';'mode radius','mean radius';'','';'median radius','mean radius'};
    error('makelogndistro:order', ...
        'For flag=%d, qt2 (%s) must be greater than qt1 (%s); got qt1=%g, qt2=%g.', ...
        flag,names{flag,2},names{flag,1},qt1,qt2);
end

rmin2=log10(rmin);
rmax2=log10(rmax);
r=logspace(rmin2,rmax2,nr);
dr=zeros(1,nr); % bin widths for the numerical integrals
dr(2:nr)=r(2:nr)-r(1:nr-1);
dr(1)=dr(2);

if flag==1
    rmod=qt1;
    rmed=qt2;
    mu=log(rmed);
    sigma=sqrt(log(rmed/rmod));
    stdev=sqrt( (exp(sigma^2)-1) * exp(2*mu+sigma^2) );
    rbar=rmed*exp(0.5*sigma^2);
elseif flag==2
    rmod=qt1;
    rbar=qt2;
    sigma=sqrt((2/3)*log(rbar/rmod));
    rmed=rmod*exp(sigma^2);
    mu=log(rmed);
    stdev=sqrt( (exp(sigma^2)-1) * exp(2*mu+sigma^2) );
elseif flag==3
    rmod=qt1;
    stdev=qt2;
    % solve rmod*(1+stdev^2/rbar^2)^(3/2)=rbar, i.e. rmod=rbar*exp(-3*sigma^2/2)
    % (the root is unique and lies in this bracket)
    ff=@(x) rmod*(1+stdev^2/x^2)^(1.5)-x;
    rbar=fzero(ff,[rmod rmod*(1+stdev^2/rmod^2)^(1.5)]);
    sigma=sqrt(log(1+(stdev/rbar)^2));
    mu=log(rbar^2/sqrt(stdev^2+rbar^2));
    rmed=exp(mu);
elseif flag==4
    rmed=qt1;
    rbar=qt2;
    mu=log(rmed);
    sigma=sqrt(2*log(rbar/rmed));
    rmod=rmed*exp(-sigma^2);
    stdev=sqrt( (exp(sigma^2)-1) * exp(2*mu+sigma^2) );
elseif flag==5
    rmed=qt1;
    stdev=qt2;
    mu=log(rmed);
    rbar=sqrt(rmed^2+rmed*sqrt(rmed^2+4*stdev^2))/sqrt(2);
    sigma=sqrt(log(1+(stdev/rbar)^2));
    rmod=exp(mu-sigma^2);
elseif flag==6
    rbar=qt1;
    stdev=qt2;
    sigma=sqrt(log(1+(stdev/rbar)^2));
    mu=log(rbar^2/sqrt(stdev^2+rbar^2));
    rmed=exp(mu);
    rmod=exp(mu-sigma^2);
elseif flag==7
    rmed=qt1;
    sigma=log(qt2);
    mu=log(rmed);
    rbar=exp(mu+0.5*sigma^2);
    rmod=exp(mu-sigma^2);
    stdev=sqrt( (exp(sigma^2)-1) * exp(2*mu+sigma^2) );
end
n=exp(-(log(r)-mu).^2/(2*sigma^2))./(r*sigma*sqrt(2*pi)); % lognormal pdf
reff=sum(r.^3.*n.*dr)/sum(r.^2.*n.*dr);
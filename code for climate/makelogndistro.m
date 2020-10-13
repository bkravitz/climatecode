function [r,n,rmod,rmed,rbar,stdev,reff]=makelogndistro(qt1,qt2,flag,rmin,rmax,nr)
% [r,n,rmod,rmed,rbar,stdev,reff]=makelogndistro(qt1,qt2,flag,rmin,rmax,nr)
%
% Ben Kravitz (ben.kravitz.work@gmail.com)  12 October 2011
%
% generates a lognormal distribution with inputs:
% - qt1 and qt2 are any two quantities from the list:
%     mode radius, median radius, mean radius, standard deviation
%     (or the combination of geometric mean and standard deviation)
%     of the lognormal distribution
% - flag specifies what qt1 and qt2 are:
%     1 - qt1=mode radius,    qt2=median radius
%     2 - qt1=mode radius,    qt2=mean radius
%     3 - qt1=mode radius,    qt2=standard deviation
%     4 - qt1=median radius,  qt2=mean radius
%     5 - qt1=median radius,  qt2=standard deviation
%     6 - qt1=mean radius,    qt2=standard deviation
%     7 - qt1=geometric mean, qt2=geometric standard deviation
%     flag=7 are commonly used quantities in reporting size distributions
% - rmin and rmax (limits of calculation of the distribution)
%     I recommend 1e-4 and 1e4, respectively
% - nr (number of radii in the interval [rmin,rmax])
%     I recommend 1e6
%
% outputs:
% - r (list of radii)
% - n (size distribution)
% - rmod (mode radius)
% - rmed (median radius, also geometric mean)
% - rbar (mean radius)
% - stdev (standard deviation)
% - reff (aerosol effective radius)

rmin2=log10(rmin);
rmax2=log10(rmax);
r=logspace(rmin2,rmax2,nr);
dr=zeros(1,nr);
dr(2:nr)=r(2:nr)-r(1:nr-1);
dr(1)=dr(2);

if flag==1;
    rmod=qt1;
    rmed=qt2;
    mu=log(rmed);
    sigma=sqrt(log(rmed/rmod));
    stdev=sqrt( (exp(sigma^2)-1) * exp(2*mu+sigma^2) );
    rbar=rmed*exp(0.5*sigma^2);
elseif flag==2;
    rmod=qt1;
    rbar=qt2;
    sigma=sqrt((2/3)*log(rbar/rmod));
    rmed=rmod*exp(sigma^2);
    mu=log(rmed);
    stdev=sqrt( (exp(sigma^2)-1) * exp(2*mu+sigma^2) );
elseif flag==3;
    rmod=qt1;
    stdev=qt2;
    ff=(1+stdev^2./r.^2).^(1.5)*rmod-r;
    I=find(ff==0);
    if isempty(I)==1;
        J=find(ff<0);
        K=find(ff>0);
        J1=intersect(J+1,K);
        K1=intersect(K+1,J);
        L=union(J1,K1);
    else
        L=I;
    end
    rbar=(r(L));
    sigma=sqrt(log(1+(stdev/rbar)^2));
    mu=log(rbar^2/sqrt(stdev^2+rbar^2));
    rmed=exp(mu);
elseif flag==4;
    rmed=qt1;
    rbar=qt2;
    mu=log(rmed);
    sigma=sqrt(2*log(rbar/rmed));
    rmod=rmed*exp(-sigma^2);
    stdev=sqrt( (exp(sigma^2)-1) * exp(2*mu+sigma^2) );
elseif flag==5;
    rmed=qt1;
    stdev=qt2;
    mu=log(rmed);
    rbar=sqrt(rmed^2+rmed*sqrt(rmed^2+4*stdev^2))/sqrt(2);
    sigma=sqrt(log(1+(stdev/rbar)^2));
    rmod=exp(mu-sigma^2);
elseif flag==6;
    rbar=qt1;
    stdev=qt2;
    sigma=sqrt(log(1+(stdev/rbar)^2));
    mu=log(rbar^2/sqrt(stdev^2+rbar^2));
    rmed=exp(mu);
    rmod=exp(mu-sigma^2);
elseif flag==7;
    rmed=qt1;
    sigma=log(qt2);
    mu=log(rmed);
    rbar=exp(mu+0.5*sigma^2);
    rmod=exp(mu-sigma^2);
    stdev=sqrt( (exp(sigma^2)-1) * exp(2*mu+sigma^2) );
end
n=lognpdf(r,mu,sigma);
nf=sum(n.*dr);
reff=sum(r.^3.*n.*dr)/sum(r.^2.*n.*dr);
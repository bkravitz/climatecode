function solarFout=insolation(ecc,eps,oma,lat,day)
% solarFout=insolation(ecc,eps,oma,lat,day)
%
% Ben Kravitz (ben.kravitz.work@gmail.com)  12 October 2011
% vectorized by Ben on 12 December 2011
% vectorization fixed by Ben on 16 February 2012
%
% This function calculates insolation received at a given
% latitude and day for specified Milankovitch parameters.
% ecc = eccentricity (range of values for Earth: 0.005-0.058)
% eps = obliquity (in radians) (22.1 to 24.5 degrees)
% oma = solar longitude of perihelion measured from the moving
%       vernal equinox (ranges in a complete circle from 0 to 2pi)
% lat = latitude (in radians)
% day = Julian day (measured from January 1)
%
% January 1 is specified to always be the winter solstice.
%
% You can modify the code directly if you want to change the
% length of the year or the solar constant.

%% References
% Berger, Andre L. (1978), Long-term variations of daily
% insolation and quaternary climatic changes, J. Atmos. Sci.,
% 35, 2362-2367,
% doi:10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2.
%
% Peixoto, Jose Pinto and Abraham H. Oort (1992),
% Physics of Climate, 520 pages.

%% parameter list
lyear=360.; % length of year (in days)
S0=1365.; % solar constant (in W/m2)

%% vectorization
[A B]=size(ecc);
if B>A;
    ecc=ecc';
end
[A,~]=size(ecc);
[C D]=size(eps);
if D>C;
    eps=eps';
end
[C,~]=size(eps);
[E F]=size(oma);
if F>E;
    oma=oma';
end
[E,~]=size(oma);
[G H]=size(lat);
if H>G;
    lat=lat';
end
[G,~]=size(lat);
[I J]=size(day);
if J>I;
    day=day';
end
[I,~]=size(day);
ecc2=zeros(A*C*E*G*I,1);
eps2=zeros(A*C*E*G*I,1);
oma2=zeros(A*C*E*G*I,1);
lat2=zeros(A*C*E*G*I,1);
for a1=1:A;
ecc2((a1-1)*C*E*G*I+1:a1*C*E*G*I)=ecc(a1);
for a2=1:C;
eps2((a1-1)*C*E*G*I+(a2-1)*E*G*I+1:a2*E*G*I+(a1-1)*C*E*G*I)=eps(a2);
for a3=1:E;
oma2((a1-1)*C*E*G*I+(a2-1)*E*G*I+(a3-1)*G*I+1:(a1-1)*C*E*G*I+(a2-1)*E*G*I+a3*G*I)=oma(a3);
for a4=1:G;
lat2((a1-1)*C*E*G*I+(a2-1)*E*G*I+(a3-1)*G*I+(a4-1)*I+1:(a1-1)*C*E*G*I+(a2-1)*E*G*I+(a3-1)*G*I+a4*I)=lat(a4);
end
end
end
end
ecc=ecc2;
eps=eps2;
oma=oma2;
lat=lat2;
clear ecc2 eps2 oma2 lat2;
day=repmat(day,[A*C*E*G 1]);

%% convenient quantities to have
ecc2=ecc.^2;
ecc3=ecc.^3;
beta=(1-ecc2).^(0.5);
OMEGA=2*pi/lyear; % angular speed of the earth
                  % also mean longitude angle subtended in one day

%% turning Jan 1 into winter solstice
a=-(pi/2);
theta=a-oma;
lmd0term1=(0.5*ecc+0.125*ecc3).*(1+beta).*sin(theta);
lmd0term2=(0.25*ecc2).*(0.5+beta).*sin(2*theta);
lmd0term3=(0.125*ecc3).*((1/3)+beta).*sin(3*theta);
lmd0=a-2*(lmd0term1-lmd0term2+lmd0term3);
day=day-lmd0/OMEGA;

%% tau=day at perihelion
tau=oma/OMEGA;

%%  This section from Berger 1978 Section 3
% l = solar longitude (angular distance along Earth's orbit)
% l==0 on the vernal equinox
% lm = mean solar longitude
% lm0 = mean solar longitude at l=0
lm0term1=(0.5*ecc+0.125*ecc3).*(1+beta).*sin(-oma);
lm0term2=-0.25*ecc2.*(0.5+beta).*sin(-2*oma);
lm0term3=0.125*ecc3.*(beta+(1/3)).*sin(-3*oma);
lm0=-2*(lm0term1+lm0term2+lm0term3);
dlm=(day-tau)*OMEGA;
lm=lm0+dlm;
l=lm+(2*ecc-0.25*ecc3).*sin(lm-oma)+1.25*ecc2.*sin(2*(lm-oma))+(13/12)*ecc3.*sin(3*(lm-oma));

%% This section from Peixoto & Oort and Berger 1978 Appendix
% delta = solar declination (lat where sun is directly overhead)
% h = hour angle from the local meridian
% Ho = hour angle at the horizon (sunrise/sunset)
rho=(1-ecc2)./(1+ecc.*cos(l-oma));
rho2=rho.^2;
delta=asin(sin(eps).*sin(l));
Ho=acos(-tan(lat).*tan(delta));
Ho( ( abs(lat) >= pi/2 - abs(delta) ) & ( lat.*delta > 0 ) )=pi;
Ho( ( abs(lat) >= pi/2 - abs(delta) ) & ( lat.*delta <= 0 ) )=0;
solarF = (S0/pi)*(Ho.*sin(lat).*sin(delta)+cos(lat).*cos(delta).*sin(Ho))./rho2;
solarFout=zeros(A,C,E,G,I);
for a1=1:A;
for a2=1:C;
for a3=1:E;
for a4=1:G;
for a5=1:I;
solarFout(a1,a2,a3,a4,a5)=solarF((a1-1)*C*E*G*I+(a2-1)*E*G*I+(a3-1)*G*I+(a4-1)*I+a5);
end
end
end
end
end
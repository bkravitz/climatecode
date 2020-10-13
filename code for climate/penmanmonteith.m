function PET=penmanmonteith(Ps,SH,LH,SAT,RH,U,height,lomask)

% function PET=penmanmonteith(Ps,SH,LH,SAT,RH,U,height,lomask)
% Coded by Ben Kravitz (ben.kravitz@pnnl.gov or ben.kravitz.work@gmail.com)
% Last updated 5 July 2017
%
% Calculates potential evapotranspiration (PET) using the Penman-Monteith equation
% Inputs:
%  Ps     = pressure at the reference height (Pa)
%  SH     = sensible heat flux (W/m2)
%  LH     = latent heat flux (W/m2)
%  SAT    = surface air temperature (K)
%  RH     = relative humidity at the reference height (unitless)
%  U      = wind speed at the reference height (m/s)
%  height = reference height (m)
%  lomask = Land/Ocean mask (=1 if land, =0 if ocean, between 0 and 1 if coastal);
%
% Outputs:
%  PET    = potential evapotranspiration (kg/m2/s)
%  Note that PET is usually expressed in m or mm as annually averaged potential
%  evapotranspiration.  To obtain this, you will need to divide by the density of
%  water (1000 kg/m3) and then sum over the appropriate times (number of seconds
%  in each timestamp that you have).
%
% References used in this code:
% Fu, Qiang and Song Feng (2014), Responses of terrestrial aridity to global warming,
%   J. Geophys. Res., 119, doi:10.1002/2014JD021608.
% Fu, Q., L. Lin, J. Huang, S. Feng, and A. Gettelman (2016), Changes in terrestrial
%   aridity for the period 850–2080 from the Community Earth System Model, J. Geophys.
%   Res. Atmos., 121, 2857–2873, doi:10.1002/ 2015JD024075.
% Richter, I., and S.-P. Xie (2008), Muted precipitation increase in global warming
%   simulations: A surface evaporation perspective, J. Geophys. Res., 113, D24118,
%   doi:10.1029/2008JD010561.
% Allen, R. G., L. S. Pereira, D. Raes, and M. Smith (1998), Crop evapotranspiration:
%   Guidelines for computing crop water requirements, FAO Irrigation and Drainage
%   Pap. 56, Stylus, Sterling, Va.
% Flatau, P. J., R. L. Walko, and W. R. Cotton (1992), Polynomial fits to
%     saturation vapor pressure, J. Appl. Met., 31, 1507-1513.
% Wexler, A. (1976), Vapor pressure formulation for water in range 0 to 100°C, A
%    revision, J. Res. Natl. Bur. Stand., 80A, 775-785.
% Scheff, J., and D. M. W. Frierson (2014), Scaling potential evapotranspiration
%    with greenhouse warming, J. Clim., 27, 1539–1558.
%

% Constants
cp = 1004; % specific heat of air (J/kg/K)
CHL = 4.8e-3; % bulk transfer coefficient (Allen et al., 1998) - land only
CHO = 1.5e-3; % bulk transfer coefficient over ocean (Richter and Xie, 2008; Fu and Song, 2014)
rsl = 70; % bulk stomatal resistance under well-watered conditions (s/m) - land only
g  = 9.81; % acceleration due to gravity (m/s2)
eps = 0.622; % molecular weight ratio of water to dry air
Lv = 2.453e6; % latent heat of vaporization (J/kg)
Rd = 287.0; % Universal gas constant for dry air (J/kg/K)
Rv = 461.51; % Universal gas constant for water (J/kg/K)
KtoC = 273.16; % Conversion factor from Kelvin to Celsius

% scaling wind from reference height to 2m height using equation 47 from Allen et al., 1998
% based on a logarithmic wind profile over a short grassy surface
u2=U*4.87/log(67.8*height-5.42);

% Delta=de*/dSAT (slope of the saturation vapor pressure curve) at temperature given by SAT
% Original code by Susan Solomon
% Clausius Clapeyron - log of saturation vapor pressure over water (Pa)
lnewsat=log(611.73)+(Lv/Rv)*(1./KtoC-1./SAT);
estar=exp(lnewsat);
Delta=estar.*Lv/Rv./(SAT.^2);

% scaling pressure from reference height to 2m height, assuming a scale height of 10 km
p2=Ps*exp((height-2)/10000);

% calculating air density
spec_humid=eps*estar./(Ps-(1-eps)*estar);
rho=p2./SAT./(Rv.*spec_humid+Rd.*(1-spec_humid));

% psychometric constant (Pa/K)
gamma = cp*p2/eps/Lv;

% creating relevant land/ocean maps (rs and CH should be different over land vs ocean)
I=size(SH);
rs=repmat(lomask*rsl,[1 1 I(3)]); % rs is 0 over oceans
CHtemp=lomask*CHL+(1-lomask)*CHO;
CH=repmat(CHtemp,[1 1 I(3)]);

% Using Penman-Monteith equation to compute PET
% Rn-G = SH+LH (Scheff and Frierson, 2014)
numerator=(SH+LH).*Delta+rho.*cp.*estar.*(1-RH/100).*CH.*abs(U);
denominator=Delta+gamma.*(1+rs.*CH.*abs(U));
PET=numerator./denominator/Lv;
PET(isnan(PET)==1)=0;
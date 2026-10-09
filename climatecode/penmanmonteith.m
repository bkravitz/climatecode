function PET=penmanmonteith(Ps,SH,LH,SAT,RH,U,height,lomask)

% function PET=penmanmonteith(Ps,SH,LH,SAT,RH,U,height,lomask)
% Coded by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 5 July 2017
% Bug fixes 8 October 2026 (with Claude):
%   (1) the 2 m wind u2 (Allen et al., 1998, eq. 47) was computed but never
%       used; over land the aerodynamic conductance is now CHL*u2, which is
%       the FAO-56 ra = 208/u2 that CHL = 4.8e-3 = 1/208 comes from (the
%       raw wind U was used before, e.g., 34% too large for U at 10 m).
%       Over ocean it remains CHO*U (Fu and Feng, 2014).
%   (2) air density now uses the actual (not saturation) specific humidity
%       at the 2 m pressure (effect < 0.5%).
%   (3) the land/ocean maps no longer require 3-D [lon lat time] inputs
%       (repmat(...,[1 1 I(3)]) replaced by implicit expansion), so
%       2-D fields, time series, and scalars work.
%   (4) missing data stay missing: NaN in any input (including lomask)
%       now gives NaN PET.  Previously all NaN PET was set to 0, so
%       missing data looked like zero evaporative demand.  There is no
%       other masking: ocean points (lomask=0) are computed as open-water
%       evaporation (rs=0), not excluded; mask them yourself if needed.
%   (5) documented the origin of CHL and CHO and the winds they are paired
%       with (see Transfer coefficients below); code unchanged.
%
% Calculates potential evapotranspiration (PET) using the Penman-Monteith
% equation in the bulk form of Scheff and Frierson (2014) and Fu and Feng
% (2014), with available energy Rn - G = SH + LH:
%   PET = [Delta*(SH+LH) + rho*cp*(es-ea)*ga] / [Delta + gamma*(1+rs*ga)] / Lv
% where ga = 1/ra is the aerodynamic conductance (m/s), rs the surface
% (stomatal) resistance, Delta = des/dT, and gamma the psychrometric
% constant.  Over land this is the FAO-56 grass reference (rs = 70 s/m,
% ra = 208/u2); over ocean rs = 0 (open water) and ga = CHO*U.  Coastal
% points (0 < lomask < 1) blend rs and ga linearly in lomask.
%
% Inputs (scalars or arrays of compatible sizes, e.g., [lon lat time] for
% the fields and [lon lat] for lomask):
%  Ps     = pressure at the reference height (Pa)
%  SH     = sensible heat flux (W/m2, positive upward, i.e., surface to
%           atmosphere, as in CESM SHFLX or CMIP hfss)
%  LH     = latent heat flux (W/m2, positive upward, e.g., CESM LHFLX or
%           CMIP hfls)
%  SAT    = surface air temperature (K)
%  RH     = relative humidity at the reference height (percent, 0-100)
%  U      = wind speed at the reference height (m/s)
%  height = reference height (m, > ~0.1; e.g., 10 for 10 m winds)
%  lomask = Land/Ocean mask (=1 if land, =0 if ocean, between 0 and 1 if coastal)
%
% Outputs:
%  PET    = potential evapotranspiration (kg/m2/s, equivalent to mm/s)
%  Multiply by 86400 to get mm/day.  Note that PET is usually expressed in
%  m or mm as annually averaged potential evapotranspiration.  To obtain
%  this, you will need to divide by the density of water (1000 kg/m3) and
%  then sum over the appropriate times (number of seconds in each
%  timestamp that you have).  PET is NaN wherever any input is NaN
%  (e.g., missing data); it is not set to 0.
%
% Transfer coefficients (ga = CH*wind)
%  Land: CHL = 4.8e-3 = 1/208 is FAO-56 (Allen et al., 1998, eq. 4 and
%  Box 4): ra = 208/u2 s/m for the 0.12 m grass reference with wind and
%  humidity measured at 2 m, so it is paired with the 2 m wind u2 (eq. 47).
%  Fu and Feng (2014) also use CH = 4.8e-3 with "u ... the wind at 2 m".
%  Ocean: CHO = 1.5e-3 is from Fu and Feng (2014, section 2.1), who cite
%  "Richter and Xie, 2008; I. Richter, personal communication, 2014".
%  Richter and Xie (2008, eq. 1) use the bulk aerodynamic formula
%  QE = Lv*CE*rho*W*(qs-qa) with W the CMIP3 model surface (10 m) wind,
%  but give no numerical CE (they hold CE constant and estimate its
%  change with the COARE algorithm of Fairall et al., 1996), so the value
%  itself rests on the personal communication.  1.5e-3 is a typical 10 m
%  neutral ocean exchange coefficient (e.g., Large and Pond, 1982, give
%  1.15e-3 for evaporation and 1.2-1.5e-3 is common), i.e., it is a bulk
%  aerodynamic coefficient for the 10 m wind; used here (as in Fu and
%  Feng) in Penman-Monteith form with rs = 0 and CH = CE.  Fu and Feng
%  describe u as the 2 m wind for both land and ocean, while CMIP output
%  winds are at 10 m; it is ambiguous which they paired with CHO.  This
%  code uses CHO with the input wind U at the reference height (so with
%  10 m winds it is the standard 10 m bulk formula).  Over open water the
%  10 m wind is ~15-20% larger than the 2 m wind (log profile), so the
%  ocean aerodynamic term differs by that much depending on the choice.
%
% Assumptions: Clausius-Clapeyron saturation vapor pressure with constant
% Lv; pressure at 2 m from a 10 km scale height; SAT and RH at the
% reference height are used as if they were 2 m values.
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
% Large, W. G., and S. Pond (1982), Sensible and latent heat flux
%   measurements over the ocean, J. Phys. Oceanogr., 12, 464-482.
% Fairall, C. W., E. F. Bradley, D. P. Rogers, J. B. Edson, and G. S. Young
%   (1996), Bulk parameterization of air-sea fluxes for TOGA COARE,
%   J. Geophys. Res., 101, 3747-3764.
% Scheff, J., and D. M. W. Frierson (2014), Scaling potential evapotranspiration
%    with greenhouse warming, J. Clim., 27, 1539–1558.
%

if any(height(:)<=0.1)
    error('penmanmonteith:height','height must be greater than 0.1 m (Allen et al., 1998, eq. 47).');
end
if all(RH(:)<=1.5) && any(RH(:)>0)
    warning('penmanmonteith:RHunits','RH appears to be a fraction; it should be in percent.');
end

% Constants
cp = 1004; % specific heat of air (J/kg/K)
CHL = 4.8e-3; % bulk transfer coefficient (Allen et al., 1998: 1/208) - land only, used with u2
CHO = 1.5e-3; % bulk transfer coefficient over ocean (Fu and Feng, 2014, citing Richter and Xie, 2008) - used with U
rsl = 70; % bulk stomatal resistance under well-watered conditions (s/m) - land only
eps = 0.622; % molecular weight ratio of water to dry air
Lv = 2.453e6; % latent heat of vaporization (J/kg)
Rd = 287.0; % gas constant for dry air (J/kg/K)
Rv = 461.51; % gas constant for water vapor (J/kg/K)
T0 = 273.16; % triple point of water (K), reference for Clausius-Clapeyron

% scaling wind from reference height to 2m height using equation 47 from Allen et al., 1998
% based on a logarithmic wind profile over a short grassy surface
u2=abs(U)*4.87./log(67.8*height-5.42);

% Delta=de*/dSAT (slope of the saturation vapor pressure curve) at temperature given by SAT
% Original code by Susan Solomon
% Clausius Clapeyron - log of saturation vapor pressure over water (Pa)
lnewsat=log(611.73)+(Lv/Rv)*(1/T0-1./SAT);
estar=exp(lnewsat);
Delta=estar.*Lv/Rv./(SAT.^2); % (Pa/K)
ea=estar.*RH/100; % actual vapor pressure (Pa)

% scaling pressure from reference height to 2m height, assuming a scale height of 10 km
p2=Ps.*exp((height-2)/10000);

% calculating air density (virtual temperature form)
spec_humid=eps*ea./(p2-(1-eps)*ea);
rho=p2./SAT./(Rv.*spec_humid+Rd.*(1-spec_humid));

% psychrometric constant (Pa/K)
gamma = cp*p2/eps/Lv;

% land/ocean dependence (rs and the aerodynamic conductance differ over land vs ocean)
rs=lomask*rsl; % surface resistance (s/m); 0 over oceans
ga=lomask.*CHL.*u2+(1-lomask).*CHO.*abs(U); % aerodynamic conductance 1/ra (m/s)

% Using Penman-Monteith equation to compute PET
% Rn-G = SH+LH (Scheff and Frierson, 2014)
numerator=(SH+LH).*Delta+rho.*cp.*(estar-ea).*ga;
denominator=Delta+gamma.*(1+rs.*ga);
PET=numerator./denominator/Lv; % NaN inputs (missing data) give NaN

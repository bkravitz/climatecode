function varargout=frostpoint(flag,varargin)

% [output1,[output2]]=frostpoint(flag,input1,input2,[units])
%
% Original code by Susan Solomon
% Modifications by Ben Kravitz (ben.kravitz.work@gmail.com) 16 February 2012
% Bug fixes 8 October 2026 (with Claude):
%   (1) flag=2: the saturation vapor pressure over water (Pa) was
%       combined with p (mb), so the mixing ratio over water was wrong
%       (usually negative); it is now computed in mb, with Murphy &
%       Koop (2005) eq. 10 replacing constant-Lv Clausius-Clapeyron
%       (which is ~37% high for supercooled water at 213 K);
%   (2) flag=2: output2 was the Flatau et al. (1992) polynomial for the
%       vapor PRESSURE (mb) over liquid WATER, not a mixing ratio over
%       ice.  It is now the saturation mixing ratio over ice, using the
%       inverse of the same frost point fit as flag=1 (so the two are
%       exact inverses).
%   (3) units: the frost point fit needs the water vapor partial pressure
%       e, which is p*water only for a VOLUME mixing ratio (mol/mol), but
%       the header documented a MASS mixing ratio (kg/kg).  An optional
%       units argument now selects either; the default is the documented
%       kg/kg, converted internally (see below).  For flag=1, this changes
%       the result for callers who were passing volume mixing ratios
%       without saying so: such callers should now pass 'mol/mol' (or
%       'ppmv').  flag=2 outputs use the same units (default kg/kg, as
%       before, with eps = 0.621981 instead of 0.622).
%
% This program can be used for two purposes:
% flag=1:  Calculate the frost point (output1, K) at a given pressure
%          (input1, mb) for given mixing ratios of water (input2, in
%          units; default kg/kg).
% flag=2:  Calculate the saturation mixing ratio (in units; default
%          kg/kg) over water (output1) and over ice (output2) for a given
%          pressure (input1, mb) and temperature (input2, K).
% Inputs may be scalars or arrays of compatible sizes.
%
% units (optional, case-insensitive; applies to input2 for flag=1 and to
% both outputs for flag=2):
%   'kg/kg'   (default) MASS mixing ratio w, kg water per kg dry air
%   'mol/mol' VOLUME mixing ratio (mole fraction) x, mol water per mol air
%   'ppmv'    volume mixing ratio in parts per million (x*1e6)
% Conversion with the vapor pressure e and pressure p (both in mb):
%   x = e/p,  w = eps*e/(p - e),  so  x = w/(eps + w) ~ w*M_air/M_h2o,
% where eps = M_h2o/M_air = 18.01528/28.9647 = 0.621981.  E.g., 5 ppmv
% = 5e-6 mol/mol = 3.11e-6 kg/kg.  If you have been passing volume
% mixing ratios (e.g., ppmv*1e-6) to flag=1, pass 'mol/mol'.
%
% Method
% Frost point (Solomon's fit): with e = p*x (mb),
%   Tf = (1.1313*ln(e) + 6182)/(24.45 - ln(e)),
% equivalently ln(e_ice) = (24.45*T - 6182)/(T + 1.1313).  This agrees with
% Murphy & Koop (2005) eq. 7 to within 0.05 K for 170-253 K, e.g.,
% 5 ppmv at 100 mb -> 192.6 K (Murphy & Koop: 192.5 K).
%
% Examples
%   Tf=frostpoint(1,100,5,'ppmv');          % 192.6 K
%   Tf=frostpoint(1,100,3.11e-6);           % same humidity in kg/kg
%   [wsw,wsi]=frostpoint(2,100,192.6,'ppmv'); % wsi = 5 ppmv
%
% References
% Murphy, D. M., and T. Koop (2005), Review of the vapour pressures of ice
%    and supercooled water for atmospheric applications, Q. J. R.
%    Meteorol. Soc., 131, 1539-1565, doi:10.1256/qj.04.94.
% (Flatau et al., 1992, and Wexler, 1976, 1977, were cited for the
%  polynomial previously used in flag=2; it is no longer used.)

if ~isscalar(flag) || ~any(flag==[1 2])
    error('frostpoint:flag','flag must be 1 or 2.');
end
if numel(varargin)~=2 && numel(varargin)~=3
    error('frostpoint:nargin','Usage: frostpoint(flag,input1,input2,[units]).');
end
units='kg/kg'; % default: mass mixing ratio, as originally documented
if numel(varargin)==3
    units=varargin{3};
end
if ~((ischar(units) && isrow(units)) || (isstring(units) && isscalar(units))) || ...
        ~any(strcmpi(units,{'kg/kg','mol/mol','ppmv'}))
    error('frostpoint:units','units must be ''kg/kg'' (default), ''mol/mol'', or ''ppmv''.');
end
units=lower(char(units));
Mh2o=18.01528; % molar mass of water (g/mol)
Mair=28.9647; % molar mass of dry air (g/mol)
epsm=Mh2o/Mair; % = 0.621981

if flag==1
    p=varargin{1}; % pressure (mb)
    water=varargin{2}; % mixing ratio of water (units)
    switch units
        case 'kg/kg'
            e=p.*water./(epsm+water); % water vapor partial pressure (mb)
        case 'mol/mol'
            e=p.*water;
        case 'ppmv'
            e=p.*water*1e-6;
    end
    lne=log(e); % ln of water vapor partial pressure (mb)
    frost=(1.1313*lne+6182)./(24.45-lne); % frost point (K)
    varargout={frost};
else
    p=varargin{1}; % pressure (mb)
    T=varargin{2}; % temperature (K)

    % Saturation vapor pressure over liquid (supercooled) water (mb)
    % Murphy & Koop (2005) eq. 10, valid 123-332 K
    ewsat=exp(54.842763-6763.22./T-4.210*log(T)+0.000367*T+ ...
        tanh(0.0415*(T-218.8)).*(53.878-1331.22./T-9.44523*log(T)+0.014025*T))/100;
    % Saturation mixing ratio over water (units)
    wsw=e2mix(ewsat,p,units,epsm);

    % Saturation vapor pressure over ice (mb): inverse of the flag=1 fit
    eisat=exp((24.45*T-6182)./(T+1.1313));
    % Saturation mixing ratio over ice (units)
    wsi=e2mix(eisat,p,units,epsm);

    varargout={wsw,wsi};
end

% If you want the values to write out to a file, feel free to
% uncomment these lines.
%
% fid=fopen('frost.out','w');
% fprintf(fid,['Frost point for ' num2str(p) ' mb\n']);
% fprintf(fid,'   Water    Frost Point\n');
% for k=1:length(water)
%     fprintf(fid,'%1.4e   %3.4f\n',water(k),frost(k));
% end
% fprintf(fid,'\n\n');
% fprintf(fid,'Saturation Mixing Ratio\n');
% fprintf(fid,'Temperature  |  Over H2O  |  Over Ice\n');
% for k=1:length(T)
%     fprintf(fid,'  %3.4f     %1.4e   %1.4e\n',T(k),wsw(k),wsi(k));
% end
% fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w=e2mix(e,p,units,epsm)
% Converts vapor pressure e (mb) at pressure p (mb) to a mixing ratio
switch units
    case 'kg/kg'
        w=epsm*e./(p-e); % mass mixing ratio
    case 'mol/mol'
        w=e./p; % volume mixing ratio (mole fraction)
    case 'ppmv'
        w=e./p*1e6;
end

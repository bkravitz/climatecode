function out=seawater_unesco(quantity,a,b,c,d,opts)
% out=seawater_unesco(quantity,...)  EOS-80 / UNESCO 1983 seawater properties
%
% rho   = seawater_unesco("density",S,T,p)
% cp    = seawater_unesco("cp",S,T,p)
% Gamma = seawater_unesco("adiabatic_lapse_rate",S,T,p)
% theta = seawater_unesco("potential_temperature",S,T,p)
% theta = seawater_unesco("potential_temperature",S,T,p,pr)
% T     = seawater_unesco("insitu_temperature",S,theta,p)
% T     = seawater_unesco("insitu_temperature",S,theta,p,pr)
% p     = seawater_unesco("pressure",z,lat)
% z     = seawater_unesco("depth",p,lat)
% ...   = seawater_unesco(...,TemperatureScale="ITS90")
%
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% originally eqnstate.m, ocncp.m, ocnpressure.m, ocntemp.m (12 October 2011)
% rewritten and consolidated 8 October 2026 (with Claude)
%
% Bugs fixed relative to the 2011 routines:
%   eqnstate:    pressure converted with press*1e-4 instead of press/10 (bar);
%                Bw coefficient k0 was 8.59035e-5 (UNESCO: 8.50935e-5)
%   ocncp:       C2 typed as 01412855 (should be 0.1412855); C2*T^2*C3*T^3
%                (should be C2*T^2+C3*T^3); e0 was -0.12331e-4 (should be
%                -1.2331e-4); pressure converted with 1e-4 instead of bar
%   ocnpressure: root selection compared pressures (Pa) with 1e4*z;
%                latitude was in radians (now DEGREES)
%   ocntemp:     pressure difference scaled by 1e-4 although the lapse rate
%                is per dbar, wrong sign, and a single Euler step (now the
%                UNESCO 4-step Runge-Kutta)
%   all:         scalar ^ and * (not vectorized)
%
% Inputs (all numeric arrays; sizes must be compatible for implicit
% expansion, e.g. a column of pressures against a row of temperatures):
%   S     = practical salinity (PSS-78)
%   T     = in-situ temperature (deg C)
%   theta = potential temperature (deg C) referenced to pr
%   p     = sea pressure (dbar; 0 at the sea surface)
%   pr    = reference pressure (dbar); default 0
%   z     = depth (m, positive down)
%   lat   = latitude (DEGREES; ocnpressure.m took radians)
% Temperatures are on the IPTS-68 scale, as in UNESCO 1983.  For ITS-90
% data either convert yourself (T68 = 1.00024*T90) or pass
% TemperatureScale="ITS90": inputs are then converted to IPTS-68 and
% temperature outputs (theta, T) converted back to ITS-90.
%
% Outputs:
%   "density"               rho (kg/m3), EOS-80 with the secant bulk
%                           modulus (pressure converted to bar internally)
%   "cp"                    specific heat at constant pressure (J/(kg K)),
%                           Millero et al. (1973) as given in UNESCO 1983
%   "adiabatic_lapse_rate"  adiabatic temperature gradient (K/dbar), Bryden (1973)
%   "potential_temperature" theta (deg C) of a parcel moved adiabatically from
%                           p to pr, Fofonoff (1977) 4-step Runge-Kutta
%   "insitu_temperature"    in-situ T (deg C) at p of a parcel with potential
%                           temperature theta referenced to pr.  Start from
%                           the UNESCO shortcut THETA(S,theta,pr,p) (error
%                           <1e-6 deg C to 4000 dbar, ~1e-4 at 10000 dbar),
%                           then iterate
%                           T <- T + (theta - THETA(S,T,p,pr)) until the
%                           round trip agrees to 1e-10 deg C.
%   "pressure"              sea pressure (dbar) from depth: Saunders (1981)
%                           closed form as first guess, then Newton
%                           iteration on UNESCO DEPTH, so it exactly
%                           inverts "depth" (to 1e-9 m)
%   "depth"                 depth (m) from pressure, UNESCO 1983 DEPTH
%                           (standard ocean S=35, T=0 deg C)
%
% Examples:
%   rho=seawater_unesco("density",35,25,10000)      % 1062.53817 kg/m3
%   th=seawater_unesco("potential_temperature",40,40,10000)  % 36.89073 C
%   T=seawater_unesco("insitu_temperature",35,2,(0:500:4000)')
%   p=seawater_unesco("pressure",(0:100:5000)',[-60 0 60])
%
% NOTE: EOS-80 / UNESCO 1983 has been superseded by TEOS-10 (IOC, SCOR and
% IAPSO, 2010: The international thermodynamic equation of seawater - 2010).
% For new work use the Gibbs SeaWater (GSW) Oceanographic Toolbox for MATLAB
% (https://www.teos-10.org/software.htm), e.g. gsw_rho, gsw_cp_t_exact,
% gsw_p_from_z, gsw_t_from_pt0.  This file is kept for consistency with
% older EOS-80-based analyses.

%% References
% Fofonoff, N. P. and R. C. Millard Jr. (1983), Algorithms for computation
% of fundamental properties of seawater, Unesco technical papers in marine
% science, 44, 53 pp.
%
% Millero, F. J., G. Perron and J. E. Desnoyers (1973), Heat capacity of
% seawater solutions from 5 to 35 C and 0.5 to 22 permil chlorinity,
% J. Geophys. Res., 78(21), 4499-4507.
%
% Bryden, H. L. (1973), New polynomials for thermal expansion, adiabatic
% temperature gradient and potential temperature of sea water, Deep-Sea
% Res., 20, 401-408.
%
% Fofonoff, N. P. (1977), Computation of potential temperature of seawater
% for an arbitrary reference pressure, Deep-Sea Res., 24, 489-491.
%
% Saunders, P. M. (1981), Practical conversion of pressure to depth,
% J. Phys. Oceanogr., 11(4), 573-574.

arguments
    quantity (1,1) string {mustBeMember(quantity,["density","cp","pressure","depth", ...
        "adiabatic_lapse_rate","potential_temperature","insitu_temperature"])}
    a {mustBeNumeric,mustBeReal}
    b {mustBeNumeric,mustBeReal}
    c {mustBeNumeric,mustBeReal} = []
    d {mustBeNumeric,mustBeReal} = []
    opts.TemperatureScale (1,1) string {mustBeMember(opts.TemperatureScale,["IPTS68","ITS90"])} = "IPTS68"
end

if ismember(quantity,["pressure","depth"])
    if ~isempty(c) || ~isempty(d)
        error('seawater_unesco:tooManyInputs','"%s" takes two inputs (z or p, and lat).',quantity);
    end
    a=double(a);
    lat=double(b);
    if quantity=="depth"
        out=unesco_depth(a,lat);
    else
        out=unesco_pressure(a,lat);
    end
    return
end

% S,T,p quantities
if isempty(c)
    error('seawater_unesco:missingPressure','"%s" requires inputs S, T, p.',quantity);
end
if isempty(d)
    d=0; % default reference pressure (dbar)
elseif ~ismember(quantity,["potential_temperature","insitu_temperature"])
    error('seawater_unesco:tooManyInputs','"%s" takes three inputs (S, T, p).',quantity);
end
S=double(a);
T=double(b);
p=double(c);
pr=double(d);
its90=opts.TemperatureScale=="ITS90";
if its90
    T=1.00024*T; % ITS-90 -> IPTS-68
end

switch quantity
    case "density"
        out=unesco_density(S,T,p);
    case "cp"
        out=unesco_cp(S,T,p);
    case "adiabatic_lapse_rate"
        out=unesco_atg(S,T,p);
    case "potential_temperature"
        out=unesco_theta(S,T,p,pr);
    case "insitu_temperature"
        % T is the potential temperature here; solve THETA(S,t,p,pr)=T for t
        th=T;
        t=unesco_theta(S,th,pr,p); % UNESCO shortcut: swap the pressures
        for k=1:20
            dt=th-unesco_theta(S,t,p,pr); % d(theta)/dt ~ 1, so fixed point converges fast
            t=t+dt;
            if all(abs(dt(:))<1e-10 | isnan(dt(:)))
                break
            end
        end
        out=t;
end

if its90 && ismember(quantity,["potential_temperature","insitu_temperature"])
    out=out/1.00024; % IPTS-68 -> ITS-90
end

end

%% local functions

function rho=unesco_density(S,T,p)
% EOS-80 density (kg/m3); S PSS-78, T IPTS-68 (deg C), p sea pressure (dbar)
P=p/10; % dbar -> bar
SR=sqrt(abs(S));
% pure water (SMOW) density
rhow=((((6.536332e-9*T-1.120083e-6).*T+1.001685e-4).*T-9.095290e-3).*T+6.793952e-2).*T+999.842594;
% one-atmosphere density
A0=(((5.3875e-9*T-8.2467e-7).*T+7.6438e-5).*T-4.0899e-3).*T+0.824493;
B0=(-1.6546e-6*T+1.0227e-4).*T-5.72466e-3;
rho0=rhow+(A0+B0.*SR+4.8314e-4*S).*S;
% secant bulk modulus K(S,T,P) (bar)
Kw=(((-5.155288e-5*T+1.360477e-2).*T-2.327105).*T+148.4206).*T+19652.21;
K0=Kw+((((-6.1670e-5*T+1.09987e-2).*T-0.603459).*T+54.6746) ...
    +((-5.3009e-4*T+1.6483e-2).*T+7.944e-2).*SR).*S;
Aw=((-5.77905e-7*T+1.16092e-4).*T+1.43713e-3).*T+3.239908;
A=Aw+((-1.6078e-6*T-1.0981e-5).*T+2.2838e-3+1.91075e-4*SR).*S;
Bw=(5.2787e-8*T-6.12293e-6).*T+8.50935e-5; % k0=8.50935e-5 (UNESCO 1983)
B=Bw+((9.1697e-10*T+2.0816e-8).*T-9.9348e-7).*S;
K=K0+(A+B.*P).*P;
rho=rho0./(1-P./K);
end

function cp=unesco_cp(S,T,p)
% specific heat (J/(kg K)), UNESCO 1983 CPSW; p in dbar
P=p/10; % dbar -> bar
SR=sqrt(abs(S));
% cp at one atmosphere
A=(-1.38385e-3*T+0.1072763).*T-7.643575;
B=(5.148e-5*T-4.07718e-3).*T+0.1770383;
C=(((2.093236e-5*T-2.654387e-3).*T+0.1412855).*T-3.720283).*T+4217.4;
cp0=(B.*SR+A).*S+C;
% pressure and temperature terms for S=0
A=(((1.7168e-8*T+2.0357e-6).*T-3.13885e-4).*T+1.45747e-2).*T-0.49592;
B=(((2.2956e-11*T-4.0027e-9).*T+2.87533e-7).*T-1.08645e-5).*T+2.4931e-4;
C=((6.136e-13*T-6.5637e-11).*T+2.6380e-9).*T-5.422e-8;
cp1=((C.*P+B).*P+A).*P;
% pressure and temperature terms for S>0
A=(((-2.9179e-10*T+2.5941e-8).*T+9.802e-7).*T-1.28315e-4).*T+4.9247e-3;
B=(3.122e-8*T-1.517e-6).*T-1.2331e-4;
A=(A+B.*SR).*S;
B=((1.8448e-11*T-2.3905e-9).*T+1.17054e-7).*T-2.9558e-6;
B=(B+9.971e-8*SR).*S;
C=(3.513e-13*T-1.7682e-11).*T+5.540e-10;
C=(C-1.4300e-12*T.*SR).*S;
cp2=((C.*P+B).*P+A).*P;
cp=cp0+cp1+cp2;
end

function G=unesco_atg(S,T,p)
% adiabatic lapse rate (K/dbar), UNESCO 1983 ATG (Bryden 1973); p in dbar
DS=S-35;
G=(((-2.1687e-16*T+1.8676e-14).*T-4.6206e-13).*p ...
    +((2.7759e-12*T-1.1351e-10).*DS+((-5.4481e-14*T+8.733e-12).*T-6.7795e-10).*T+1.8741e-8)).*p ...
    +(-4.2393e-8*T+1.8932e-6).*DS ...
    +((6.6228e-10*T-6.836e-8).*T+8.5258e-6).*T+3.5803e-5;
end

function th=unesco_theta(S,T0,p0,pr)
% potential temperature (deg C) at pr of a parcel at (S,T0,p0), UNESCO 1983
% THETA: Fofonoff (1977) 4th-order Runge-Kutta (Gill) in one step from p0 to pr
h=pr-p0;
xk=h.*unesco_atg(S,T0,p0);
T=T0+0.5*xk;
q=xk;
p=p0+0.5*h;
xk=h.*unesco_atg(S,T,p);
T=T+0.29289322*(xk-q);
q=0.58578644*xk+0.121320344*q;
xk=h.*unesco_atg(S,T,p);
T=T+1.707106781*(xk-q);
q=3.414213562*xk-4.121320344*q;
p=p+0.5*h;
xk=h.*unesco_atg(S,T,p);
th=T+(xk-2*q)/6;
end

function z=unesco_depth(p,lat)
% depth (m) from sea pressure p (dbar) and latitude (degrees), UNESCO 1983 DEPTH
x=sind(lat).^2;
gr=9.780318*(1+(5.2788e-3+2.36e-5*x).*x)+1.092e-6*p; % gravity incl. pressure term
z=((((-1.82e-15*p+2.279e-10).*p-2.2512e-5).*p+9.72659).*p)./gr;
end

function p=unesco_pressure(z,lat)
% sea pressure (dbar) from depth z (m) and latitude (degrees): Saunders (1981)
% closed form, refined by Newton iteration so that it inverts unesco_depth
c1=(5.92+5.25*sind(lat).^2)*1e-3;
p=((1-c1)-sqrt((1-c1).^2-8.84e-6*z))/4.42e-6; % Saunders (1981), within ~0.1%
x=sind(lat).^2;
g0=9.780318*(1+(5.2788e-3+2.36e-5*x).*x);
for k=1:20
    gr=g0+1.092e-6*p;
    zm=(((-1.82e-15*p+2.279e-10).*p-2.2512e-5).*p+9.72659).*p;
    dzm=((-4*1.82e-15*p+3*2.279e-10).*p-2*2.2512e-5).*p+9.72659;
    f=zm./gr-z;
    dp=f./((dzm.*gr-zm*1.092e-6)./gr.^2);
    p=p-dp;
    if all(abs(f(:))<1e-9 | isnan(f(:)))
        break
    end
end
end

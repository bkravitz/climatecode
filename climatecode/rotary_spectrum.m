function [f,Scw,Sccw,Stot,rc] = rotary_spectrum(u,v,dt,opts)
% [f,Scw,Sccw,Stot,rc] = rotary_spectrum(u,v,dt)
% [f,Scw,Sccw,Stot,rc] = rotary_spectrum(u,v,dt,Detrend=..,Window=..)
%
% Rotary spectra of a 2-D vector time series (e.g. currents or winds):
% splits the variance of w = u + i*v into clockwise and counterclockwise
% rotating components as a function of frequency (Gonella 1972; Emery and
% Thomson, Data Analysis Methods in Physical Oceanography, Section 5.8).
%
% Inputs:
%   u, v  = east and north components, vectors of the same length N,
%           uniformly sampled, no NaNs
%   dt    = sampling interval (optional, default 1); frequencies are in
%           cycles per unit of dt (e.g. dt=1/24 for hourly data in days
%           gives cycles per day)
%   Detrend = "mean" (default), "linear" or "none", applied to u and v
%   Window  = "none" (default, raw periodogram) or "hann"; with a window
%           the spectra are normalized by the window energy so the
%           variance is preserved on average
%
% Outputs (column vectors, f = 0 to the Nyquist frequency):
%   f     = frequency (cycles per unit time), (0:floor(N/2))'/(N*dt)
%   Scw   = clockwise power spectral density
%   Sccw  = counterclockwise power spectral density
%           (units of u^2 per unit frequency, e.g. (m/s)^2/cpd)
%   Stot  = Scw + Sccw, the total (u plus v) spectrum; sum(Stot)/(N*dt)
%           equals the variance of u plus the variance of v (Parseval;
%           exact with Window="none")
%   rc    = rotary coefficient (Sccw - Scw)./(Sccw + Scw), from -1 (purely
%           clockwise) to +1 (purely counterclockwise)
%
% Method: with MATLAB's fft convention, a counterclockwise rotating vector
% exp(+i*2*pi*f*t) appears at positive frequency in fft(u+i*v) and a
% clockwise one at negative frequency, so Sccw(f) = P(+f) and Scw(f) =
% P(-f), with P the two-sided periodogram of w.  At f=0 (and at the Nyquist
% frequency for even N) the two signs are the same bin, so its power is
% split equally between Scw and Sccw.  With x east and y north,
% counterclockwise is cyclonic in the Northern Hemisphere.
%
% These are raw (unsmoothed) estimates with 2 degrees of freedom each; to
% reduce the noise, band-average them, e.g. [Scw_s,df] = tsmooth(Scw,kw).
%
% Written 9 October 2026 by Ben Kravitz (bkravitz@iu.edu or
% ben.kravitz.work@gmail.com) with Claude, replacing the course scripts
% rotary.m, rotary2.m and rotary3.m (after Bob Chant).  rotary2.m did not
% separate the two rotation senses and plotted against a misaligned
% frequency axis; rotary3.m gave correct labels only for column input.
%
% References:
% Gonella, J. (1972), A rotary-component method for analysing meteorological
% and oceanographic vector time series, Deep-Sea Res., 19, 833-846.
% Emery, W. J. and R. E. Thomson (2001), Data Analysis Methods in Physical
% Oceanography, 2nd ed., Elsevier.

arguments
    u (:,1) double {mustBeFinite}
    v (:,1) double {mustBeFinite}
    dt (1,1) double {mustBePositive} = 1
    opts.Detrend (1,1) string {mustBeMember(opts.Detrend,["mean","linear","none"])} = "mean"
    opts.Window (1,1) string {mustBeMember(opts.Window,["none","hann"])} = "none"
end
if numel(u) ~= numel(v)
    error('rotary_spectrum:size','u and v must have the same length.');
end
N = numel(u);

% remove the mean or trend from each component
switch opts.Detrend
    case "mean"
        u = u - mean(u);
        v = v - mean(v);
    case "linear"
        u = detrend(u);
        v = detrend(v);
end

% taper; U = window energy so that a white spectrum keeps its level
if opts.Window == "hann"
    win = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);
else
    win = ones(N,1);
end
U = sum(win.^2);

% two-sided periodogram of w = u + iv
P = dt/U*abs(fft((u + 1i*v).*win)).^2;

% positive frequencies -> counterclockwise, negative -> clockwise
nf = floor(N/2) + 1;
k = (0:nf-1)';
f = k/(N*dt);
Sccw = P(k+1);
Scw = P(mod(-k,N)+1);
shared = (k == 0) | (2*k == N); % f=0 and (even N) Nyquist: one bin for both
Sccw(shared) = Sccw(shared)/2;
Scw(shared) = Scw(shared)/2;

Stot = Scw + Sccw;
rc = (Sccw - Scw)./Stot;
end

function [kp,ki,negkp] = pi_gains(wgc,amplitude_ratio,phase_lag,margins,opts)
% [kp,ki,negkp] = pi_gains(wgc,amplitude_ratio,phase_lag,margins)
% [kp,ki,negkp] = pi_gains(...,Delay=0.5)
%
% Proportional-integral (PI) controller gains that place the gain crossover
% at frequency wgc with specified phase margins, given the plant's
% amplitude ratio and phase lag at wgc.  This is a computation tool: it
% returns the gains for every requested margin and does not choose one.
%
% Controller:  C(s) = kp + ki/s.  The gains satisfy |C(i*wgc)*G(i*wgc)| = 1
% and the requested phase margin, from
%    alpha = margin - 90 + lag
%    ki = wgc/amplitude_ratio/sqrt(1+tan(alpha)^2),  kp = ki/wgc*tan(alpha)
%
% Inputs:
%   wgc             = gain crossover frequency (rad per unit time, e.g.
%                     2*pi/32 for a 32-year period with time in years)
%   amplitude_ratio = plant amplitude ratio |G| at wgc
%   phase_lag       = plant phase lag at wgc (degrees, positive = lag)
%   margins         = phase margins to design for (degrees), vector;
%                     default 5:5:90
%   Delay           = (name-value, default 0.5) delay in units of time
%                     steps used in the half-time-step phase correction:
%                     lag = phase_lag - Delay*wgc*180/pi.  The correction
%                     and its sign follow Kravitz et al. (2016); set
%                     Delay=0 to use phase_lag as given.
%
% Outputs (column vectors, one row per margin):
%   kp, ki  = proportional and integral gains
%   negkp   = true where kp < 0.  Those margins are flagged, not removed:
%             a negative proportional gain is a valid solution of the
%             equations but usually not a practical controller.
% With no outputs, prints a table of margin, kp, ki, kp/ki and the flag.
%
% Sign convention: the gains follow the loop sign convention of the
% feedback design in Kravitz et al. (2016), in which the sign of the plant
% (e.g. more SRM forcing gives lower temperature) is taken up by the
% feedback loop rather than by these formulas.  Treated as a plain positive
% plant G = amplitude_ratio*exp(-i*lag), the same gains give a phase margin
% exactly 180 degrees from the requested one, so use them only with that
% convention.
%
% Example (the values from the original computegains.m script):
%   pi_gains(2*pi/32, 8.4108, 216.0959, 5:5:90)
%
% Reference:
% Kravitz, B., D. G. MacMartin, H. Wang, and P. J. Rasch (2016),
% Geoengineering as a design problem, Earth Syst. Dynam., 7, 469-497,
% doi:10.5194/esd-7-469-2016.
%
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Converted from the script computegains.m to a function 9 October 2026
% (with Claude); the formulas are unchanged.

arguments
    wgc (1,1) double {mustBePositive}
    amplitude_ratio (1,1) double {mustBePositive}
    phase_lag (1,1) double {mustBeReal}
    margins (:,1) double {mustBeReal} = (5:5:90)'
    opts.Delay (1,1) double {mustBeNonnegative} = 0.5
end

lag = phase_lag - opts.Delay*wgc*180/pi; % half-time-step correction (Kravitz et al. 2016)
alpha = margins - 90 + lag;
ki = wgc/amplitude_ratio./sqrt(1+tand(alpha).^2);
kp = ki/wgc.*tand(alpha);
negkp = kp < 0;

if nargout == 0
    fprintf('%8s %12s %12s %12s  %s\n','margin','kp','ki','kp/ki','');
    for n = 1:numel(margins)
        flag = '';
        if negkp(n)
            flag = 'kp < 0';
        end
        fprintf('%8.1f %12.5g %12.5g %12.5g  %s\n',margins(n),kp(n),ki(n),kp(n)/ki(n),flag);
    end
    clear kp
end
end

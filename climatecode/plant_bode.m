function [models,gains] = plant_bode(opts)
% [models,gains] = plant_bode
% [models,gains] = plant_bode(Name=Value,...)
%
% Compares three simple plant models of a climate response that all have
% the same gain and phase at the design frequency om, designs a PI
% controller for each that gives the requested phase margin at om, and
% (optionally) draws Bode plots of the plants G and the loops GK and the
% sensitivity |1/(1+GK)|.
%
% Plant models (the ARX(1) model is specified; the others are fitted to it
% at om; all include the same pure delay exp(-D*s)):
%   ARX(1)  G1 = exp(-D*s)*beta/(1+tau1*s), with tau1 chosen so that the
%           phase of beta/(1+tau1*s) at om is Phase+PhaseAdded
%   SD      Gd = exp(-D*s)*betad/(1+sqrt(taud*s))  (semi-infinite diffusion)
%   O2      G2 = exp(-D*s)*beta2/((tau2*s)^2+2*zeta*tau2*s+1)
% The SD model lags by at most 45 degrees, so it can only be fitted when
% Phase+PhaseAdded is between -45 and 0 degrees; outside that range it is
% returned as NaN with a warning.
%
% Name-value inputs (defaults are the values of the original megabode.m):
%   om           = design (gain crossover) frequency, rad/year (0.2)
%   Phase        = ARX(1) plant phase at om, degrees, negative = lag (-27)
%   Beta         = ARX(1) zero-frequency gain (0.447)
%   D            = pure delay, years (0)
%   Zeta         = damping ratio of the second-order model (2)
%   PhaseAdded   = extra phase added to Phase before fitting, degrees (0);
%                  e.g. 0.5*om*180/pi for a half-time-step correction
%   PhaseMargin  = desired phase margin, degrees (90)
%   w            = frequencies for the curves, rad/year (0.01:0.01:10)
%   Plot         = draw the three figures (true)
%
% Outputs:
%   models  = struct with fields
%       w                   frequency vector (row)
%       G1, Gd, G2          complex frequency responses (rows)
%       tau1, taud, betad, tau2, beta2   fitted parameters
%       GAtOm               complex responses [G1 Gd G2] exactly at om
%       gainAtOm, phaseAtOm gain and phase (deg) of the plants at om
%                           (identical for the three models by construction)
%   gains   = struct with fields kp and ki, each [ARX(1) SD O2], from
%             pi_gains(om, gainAtOm, -phaseAtOm, PhaseMargin, Delay=0);
%             also the loops H1, Hd, H2 (PI) and sensitivities S1, Sd, S2
%
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Converted from the script megabode.m to a function 9 October 2026 (with
% Claude).  With the default inputs the curves and figures are the same as
% the script's.  Changes: the responses are evaluated exactly at om (the
% script searched the frequency grid, which could return zero or two
% points for other om); the delay is now applied to all three models (the
% script left it off the second-order model); the PI gains use each
% model's actual gain and phase at om, including the delay and PhaseAdded
% (the script used Phase alone, so the requested margin was missed when
% D or PhaseAdded was nonzero); the SD fit now warns when it has no
% solution.
%
% See also PI_GAINS.

arguments
    opts.om (1,1) double {mustBePositive} = 0.2
    opts.Phase (1,1) double {mustBeReal} = -27
    opts.Beta (1,1) double {mustBePositive} = 0.447
    opts.D (1,1) double {mustBeNonnegative} = 0
    opts.Zeta (1,1) double {mustBePositive} = 2
    opts.PhaseAdded (1,1) double {mustBeReal} = 0
    opts.PhaseMargin (1,1) double {mustBeReal} = 90
    opts.w (1,:) double {mustBePositive} = 0.01:0.01:10
    opts.Plot (1,1) logical = true
end
om = opts.om; D = opts.D; zeta = opts.Zeta; w = opts.w;
s = 1i*w;
som = 1i*om;
delay = @(x) exp(-D*x);

% ARX(1): phase of the rational part at om is Phase+PhaseAdded
plantphase = opts.Phase + opts.PhaseAdded;
tau1 = -tand(plantphase)/om;
arx = @(x) opts.Beta./(1+x*tau1);
G1 = delay(s).*arx(s);
gainAtOm = abs(delay(som)*arx(som));
ratphase = angle(arx(som))*180/pi; % phase of G1 at om without the delay

% SD: phase of 1/(1+sqrt(i*taud*om)) is -atan(a/(1+a)), a = sqrt(taud*om/2)
t = tand(ratphase);
if t <= -1 || t > 0
    warning('plant_bode:SDphase', ['The semi-infinite diffusion model lags ' ...
        'by at most 45 degrees and cannot reach a phase of %.1f degrees; ' ...
        'SD results are NaN.'], ratphase);
    taud = NaN;
else
    taud = (t/(1+t))^2*2/om;
end
betad = gainAtOm*abs(1+sqrt(som*taud));
sd = @(x) betad./(1+sqrt(taud*x));
Gd = delay(s).*sd(s);

% O2: smaller positive root of tan(phase)*x^2 - 2*zeta*x - tan(phase) = 0,
% x = tau2*om
tau2 = -(-2*zeta*om+sqrt(4*om^2*(zeta^2+t^2)))/(2*om^2*t);
beta2 = gainAtOm*abs((tau2*som)^2+2*zeta*tau2*som+1);
o2 = @(x) beta2./((tau2*x).^2+2*zeta*tau2*x+1);
G2 = delay(s).*o2(s);

% PI gains from each model's gain and phase at om (delay included)
Gom = [delay(som)*arx(som), delay(som)*sd(som), delay(som)*o2(som)];
phaseAtOm = angle(Gom(1))*180/pi;
kp = NaN(1,3); ki = NaN(1,3);
for n = 1:3
    if ~isfinite(Gom(n)) % SD model could not be fitted
        continue
    end
    [kp(n),ki(n)] = pi_gains(om,abs(Gom(n)),-angle(Gom(n))*180/pi, ...
        opts.PhaseMargin,Delay=0);
end

H1 = G1.*(ki(1)./s+kp(1));  Hd = Gd.*(ki(2)./s+kp(2));  H2 = G2.*(ki(3)./s+kp(3));
models = struct('w',w,'G1',G1,'Gd',Gd,'G2',G2,'tau1',tau1,'taud',taud, ...
    'betad',betad,'tau2',tau2,'beta2',beta2,'GAtOm',Gom,'gainAtOm',gainAtOm, ...
    'phaseAtOm',phaseAtOm);
gains = struct('kp',kp,'ki',ki,'H1',H1,'Hd',Hd,'H2',H2, ...
    'S1',1./(1+H1),'Sd',1./(1+Hd),'S2',1./(1+H2));

if opts.Plot
    plotFigures(models,gains,om);
end
end

%% ------------------------------------------------------------------------
function plotFigures(m,g,om)
% The three figures of the original megabode.m
w = m.w; s = 1i*w;
deg = @(x) angle(x)*180/pi;
H1a = m.G1.*g.ki(1)./s;  Hda = m.Gd.*g.ki(2)./s;  H2a = m.G2.*g.ki(3)./s;
xt = [0.01:0.01:0.09 0.1:0.1:0.9 1:10];
xtl = {'1/100','1/50','','','1/20','','','','','1/10','1/5','','','1/2', ...
    '','','','','1','2','','','5','','','','','10'};
gainOm = m.gainAtOm; phaseOm = m.phaseAtOm;

% Bode plot of G
figure;
subplot(2,1,1);
loglog(w,abs(m.G1),'b','LineWidth',2); hold on;
loglog(w,abs(m.Gd),'r','LineWidth',2);
loglog(w,abs(m.G2),'g','LineWidth',2);
legend({'ARX(1)','SD','Order 2'},'Location','SouthWest','FontSize',15,'Box','off','Color','none');
set(gca,'FontSize',15,'XTick',xt,'XTickLabel',xtl,'YTick',xt, ...
    'YTickLabel',{'0.01','','','','','','','','','0.1','','','','','','','','','1','','','','','','','','','10'});
ylabel('Magnitude(G)','FontSize',15);
xlim([0.01 1.5]); ylim([0.1 1]);
plot([0.01 10],[gainOm gainOm],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
plot([om om],[0.01 gainOm],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
xlabel('Angular Frequency (rad/year)','FontSize',15);
subplot(2,1,2);
semilogx(w,deg(m.G1),'b','LineWidth',2); hold on;
semilogx(w,deg(m.Gd),'r','LineWidth',2);
semilogx(w,deg(m.G2),'g','LineWidth',2);
xlabel('Angular Frequency (rad/year)','FontSize',15);
set(gca,'FontSize',15,'XTick',xt,'XTickLabel',xtl,'YTick',-180:30:0);
ylim([-180 0]); xlim([0.01 1.5]);
ylabel('Phase(G) (degrees)','FontSize',15);
plot([0.01 10],[phaseOm phaseOm],'--','Color',[0.5 0.5 0.5]);
plot([om om],[-180 0],'--','Color',[0.5 0.5 0.5]);

% Bode plot of GK (integral-only and PI)
figure;
subplot(2,1,1);
loglog(w,abs(H1a),'k','LineWidth',2); hold on;
loglog(w,abs(g.H1),'r','LineWidth',2);
loglog(w,abs(Hda),'g','LineWidth',2);
loglog(w,abs(g.Hd),'b','LineWidth',2);
loglog(w,abs(H2a),'c','LineWidth',2);
loglog(w,abs(g.H2),'m','LineWidth',2);
set(gca,'FontSize',15,'XTick',xt,'XTickLabel',xtl,'YTick',[0.05:0.01:0.09 0.1:0.1:0.9 1:10], ...
    'YTickLabel',{'','','','','','0.1','','','','','','','','','1','','','','','','','','','10'});
ylabel('Magnitude(GK)','FontSize',15);
xlim([0.01 0.6]); ylim([0.1 30]);
legend({'ARX(1) (I)','ARX(1) (PI)','SD (I)','SD (PI)','O2 (I)','O2 (PI)'},'Location','SouthWest','FontSize',15,'Box','off','Color','none');
xlabel('Angular Frequency (rad/year)','FontSize',15);
subplot(2,1,2);
semilogx(w,deg(H1a),'k','LineWidth',2); hold on;
semilogx(w,deg(g.H1),'r','LineWidth',2);
semilogx(w,deg(Hda),'g','LineWidth',2);
semilogx(w,deg(g.Hd),'b','LineWidth',2);
semilogx(w,deg(H2a),'c','LineWidth',2);
semilogx(w,deg(g.H2),'m','LineWidth',2);
xlabel('Angular Frequency (rad/year)','FontSize',15);
set(gca,'FontSize',15,'XTick',xt,'XTickLabel',xtl,'YTick',-210:30:0);
ylim([-185 0]); xlim([0.01 0.6]);
ylabel('Phase(GK) (degrees)','FontSize',15);
plot([om om],[-210 0],'--','Color',[0.5 0.5 0.5]);
plot([0.01 2.5],[-180 -180],'k');

% Sensitivity function
figure;
semilogx(w,abs(g.S1),'b','LineWidth',2); hold on;
semilogx(w,abs(g.Sd),'r','LineWidth',2);
semilogx(w,abs(g.S2),'g','LineWidth',2);
set(gca,'FontSize',15,'XTick',xt,'XTickLabel',xtl);
ylabel('Sensitivity','FontSize',15);
xlim([0.01 2.5]); ylim([0 1.4]);
S1om = abs(1/(1+m.GAtOm(1)*(g.ki(1)/(1i*om)+g.kp(1)))); % |S| of the ARX(1) loop at om
plot([0.01 10],[S1om S1om],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
plot([om om],[0 S1om],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
plot([0.01 10],[1 1],'-','Color','k','HandleVisibility','off');
text(1/75,1.04,'Amplification','FontSize',13);
text(1/75,0.96,'Attenuation','FontSize',13);
legend({'ARX(1)','SD','O2'},'FontSize',15,'Location','SouthEast','Box','off');
xlabel('Angular Frequency (rad/year)','FontSize',15);
end

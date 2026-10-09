function [ma,mi,ang]=pca_ellipse(un,ue)
% PCA_ELLIPSE  Principal axes (variance ellipse) of a 2-D vector time series
%function [ma,mi,ang]=pca_ellipse(un,ue)
%
% Computes the principal-axis (variance ellipse) decomposition of a pair of
% velocity components, e.g. ocean currents or winds: the eigenvalues and the
% major-axis orientation of the 2x2 covariance matrix [var(ue) cov; cov var(un)]
% (see e.g. Emery and Thomson, Data Analysis Methods in Physical
% Oceanography, section on principal axes).
%
% Inputs:
%   un = northward (y) component, vector
%   ue = eastward (x) component, vector of the same length as un
%   (any consistent units; NaNs are not removed and will propagate)
%
% Outputs:
%   ma  = variance along the major axis (units of input squared)
%   mi  = variance along the minor axis (units of input squared)
%         (use sqrt(ma), sqrt(mi) for the ellipse semi-axis lengths)
%   ang = orientation of the major axis in radians, measured counterclockwise
%         from east (the ue axis), in the range (-pi/2, pi/2]
%
% Variances use the 1/N (population) normalization.
%
% Renamed from pca.m on 8 October 2026 to avoid shadowing the Statistics
% Toolbox function pca.
% Bug fixes 8 October 2026 (with Claude): major-axis angle now computed as
%   atan2(2*cov,var(ue)-var(un))/2, which is defined when ma==var(un) (the old
%   atan(cov/(ma-var(un))) gave NaN, e.g. for uncorrelated components with
%   var(un)>var(ue), where the answer is pi/2); local variable cov renamed so
%   it no longer shadows the built-in cov.

% remove means (so the mean(.)*mean(.) terms below are ~0)
ue=ue-mean(ue);
un=un-mean(un);
u2=ue.*ue;
v2=un.*un;
uv=ue.*un;
usd=mean(u2)-mean(ue)*mean(ue);
vsd=mean(v2)-mean(un)*mean(un);
uvcov=mean(uv)-mean(ue)*mean(un);
D=(usd+vsd)^2;
E=sqrt(D-4*(usd*vsd-uvcov*uvcov)); % eigenvalues of the covariance matrix are (trace +/- E)/2
ma=(usd+vsd+E)/2;
mi=(usd+vsd-E)/2;
ang=atan2(2*uvcov,usd-vsd)/2; % same as atan(uvcov/(ma-vsd)), but defined when ma==vsd
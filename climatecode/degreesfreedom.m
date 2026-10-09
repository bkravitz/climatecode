function [dof,Tscale,xc,lags] = degreesfreedom(D,sided)
% [dof,Tscale,xc,lags] = degreesfreedom(D);
% [dof,Tscale,xc,lags] = degreesfreedom(D,sided);
% Estimate degrees of freedom and autocorrelation scale for the ensemble 
% of time series in rows of data matrix D.
%
% The method assumes the rows of D are a set of time series of 
% length N = size(D,2) at M = size(D,1) different locations.
%
% Inputs:
%   D = M x N data matrix (rows are time series, uniformly sampled)
%   sided = (optional) form of the integral autocorrelation time scale,
%           with r(k) the ensemble-mean autocorrelation at lag k, summed
%           over lags k >= 1 before the first zero crossing:
%      'one-sided' or 1 (default):  Tscale = 1 + sum r(k)
%      'two-sided' or 2:  Tscale = 1 + 2*sum r(k)  (Trenberth 1984;
%           Bretherton et al. 1999).  For an AR(1) process with lag-1
%           autocorrelation rho this tends to (1+rho)/(1-rho).
%      'ar1' or 3:  Tscale = (1+r(1))/(1-r(1)), assuming an AR(1) process
%           and using only the ensemble-mean lag-1 autocorrelation r(1)
%           (Bretherton et al. 1999, eq. 31; Zwiers and von Storch 1995).
%           N/Tscale is then the usual AR(1) effective sample size
%           N(1-r(1))/(1+r(1)).  This is more stable than the summed
%           estimates for short time series, where the sample
%           autocorrelation at larger lags (and the first zero crossing)
%           is noisy; the summed estimates are better when the series is
%           not AR(1)-like (e.g. oscillatory or long-memory).  If r(1) < 0,
%           Tscale < 1 and dof exceeds N.
%   Strings are case-insensitive.
%
% Outputs:
%   dof = degrees of freedom, round(N/Tscale) (use N/Tscale for the
%         unrounded effective sample size)
%   Tscale = autocorrelation time scale (in units of the sampling interval)
%   xc = ensemble-mean autocorrelation at lags 0,1,...,N-1 (xc(1)=1)
%   lags = the lags of xc (0,1,...,N-1)
%
% John Wilkin, April 2004
% Optional 'sided' input (one- or two-sided time scale) added 8 October 2026
% 'ar1' option added 9 October 2026 (with Claude); replaces the separate
%   effectivesamplesize.m (single-series AR(1) formula, which needed the
%   Econometrics Toolbox)
%
% Bretherton, C. S., M. Widmann, V. P. Dymnikov, J. M. Wallace, and
% I. Blade (1999), The effective number of spatial degrees of freedom of a
% time-varying field, J. Climate, 12, 1990-2009.
% Trenberth, K. E. (1984), Some effects of finite sample size and
% persistence on meteorological statistics.  Part I: Autocorrelations,
% Mon. Wea. Rev., 112, 2359-2368.
% Zwiers, F. W. and H. von Storch (1995), Taking serial correlation into
% account in tests of the mean, J. Climate, 8, 336-351.

% Check the 'sided' option
if nargin < 2 || isempty(sided)
   sided = 1;
end
names = {'one-sided','two-sided','ar1'};
if (ischar(sided) || (isstring(sided) && isscalar(sided))) && ...
      any(strcmpi(sided,names))
   method = find(strcmpi(sided,names));
elseif isnumeric(sided) && isscalar(sided) && any(sided == [1 2 3])
   method = sided;
else
   error('degreesfreedom:badSided', ['sided must be ''one-sided'' (or 1), ' ...
      '''two-sided'' (or 2) or ''ar1'' (or 3).'])
end

% Compute the covariance of each time series
%   If rows of D are all zero (or constant) the output will be NaN
xc = NaN(size(D,1),2*size(D,2)-1);
for m=1:size(D,1)
   % 'coeff' gives normalized function xc(0)=1
   [xc(m,:),lags] = xcov(D(m,:),'coeff');
end

% remove the rows with NaNs
xc(any(isnan(xc),2),:) = [];

% Take mean over all M time series
xc = mean(xc,1);

% Retain half of cross-correlation sequence that starts from lag=0
xc = xc(find(lags==0):end);
lags = lags(find(lags==0):end);

% Integrate to first zero crossing, at t=t*, and find autocorrelation scale 
% Tscale such that: Tscale = integral_0^t* xc(t') dt' 
% (This assumes time series are recorded at uniform sampling intervals.)
% Tscale will be in the normalized units of the sampling interval.
% One-sided:  Tscale = xc(0) + sum xc(k), k = 1 to t*-1
% Two-sided:  lags k = 1 to t*-1 counted for both +k and -k
% AR(1):      Tscale = (1+r(1))/(1-r(1)) from the lag-1 value r(1) = xc(2)
if method == 3
   Tscale = (1+xc(2))/(1-xc(2)); % xc(2) is lag 1
else
   tstar = find(xc<=0,1)-1;
   Tscale = mean(xc(1:tstar))*tstar;
   if method == 2
      Tscale = 2*Tscale - xc(1); % = xc(0) + 2*sum xc(k)
   end
end

% Degrees of Freedom is estimated as the number of correlation time scales 
% in length of the data record
dof = round(size(D,2)/Tscale);

function outfield=scalespace(data,sigma,N)

% function outfield=scalespace(data,sigma,N)
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 17 September 2014
% Bug fixes 8 October 2026 (with Claude): the spherical harmonic degree
%   used in the diffusion filter was off by one (degree l was damped as if
%   it were l+1), which over-smoothed every scale and damped the global mean
%   by exp(-sigma^2); the global mean is now preserved.  The output is now
%   real (round-off imaginary parts dropped).  These change the results.
% Changes 8 October 2026 (with Claude): the dependency on SPHRHARM
%   (VOICEBOX, Mike Brookes) has been removed.  The spherical harmonic
%   filter is now computed here directly (independent implementation, same
%   grid and normalization conventions), so this file is self-contained.
%   Results are identical to the SPHRHARM version to round-off.
%
% This function computes scale space methods (e.g., Marvel et al., 2014, JGR)
% using smoothing via diffusion on a sphere.  This has been shown to be a
% more accurate way of comparing 2-D climate model output by removing small
% scale features.  This is better than averaging because some features in
% models may be due to persistent biases, which would affect the averaged
% ensemble results.
%
% Inputs:  data  = The data you are attempting to smooth.  Transform it so that
%                  the matrix has size (number of lats) x (number of lons)
%                  It must be a global, regularly spaced lat-lon grid with
%                  latitudes at the centres of equal bands, i.e.
%                  -90+(k-1/2)*180/nlat for k=1:nlat (e.g. -88.75:2.5:88.75),
%                  and longitudes equally spaced around the full circle.
%                  Either latitude order (N-S or S-N) and any starting
%                  longitude work, since the filter is symmetric.
%          sigma = The spatial scale of smoothing you want to do.  Must be in
%                  radians.  For example, if you want to do 2 degree smoothing,
%                  put sigma=2*pi/180.
%          N     = The spectral truncation for the spherical harmonics computation.
%                  In theory, higher N values are more accurate but more
%                  computationally expensive.  In practice, transforming a 2-D
%                  field is not that computationally expensive, so you can pretty
%                  much put whatever you want.  I would recommend using values
%                  between 20 and 100, depending on how much fine scale structure
%                  is in your data set.  Degrees above N are removed entirely.
%                  For an exact transform the grid needs nlat >= 2N+1 and
%                  nlon >= 2N+1; a warning is given otherwise.
%
% Outputs:  outfield = The smoothed data set, same size as data (lats x lons)
%
% Method: diffusion on the unit sphere, applied in spectral space.  Each
% spherical harmonic coefficient of degree l is multiplied by
% exp(-l*(l+1)*sigma^2/2), the heat kernel for diffusion time sigma^2/2.
% The transform is done one zonal wavenumber m at a time: an FFT along
% longitude, then projection onto the fully normalized associated Legendre
% functions P_l^m(cos(colatitude)), l=m:N, using Fejer's first quadrature
% rule in latitude (exact for the band-centre grid above), damping, and
% synthesis back to the grid.

[nlat,nlon]=size(data);
if 2*N+1>nlat || 2*N+1>nlon
    warning('scalespace:coarsegrid',['Grid (%d x %d) is too coarse for an ' ...
        'exact transform at N=%d (need at least %d x %d).'],nlat,nlon,N,2*N+1,2*N+1);
end

% colatitudes of the band centres and Fejer (first rule) quadrature weights
% for integrating over mu=cos(colatitude) on [-1,1]
theta=((1:nlat)'-0.5)*pi/nlat;
j=1:floor(nlat/2);
w=(2/nlat)*(1-2*cos(2*theta*j)*(1./(4*j'.^2-1)));

% fully normalized associated Legendre functions: P(l+1,m+1,:) = P_l^m(mu),
% with the integral of P_l^m^2 over [-1,1] equal to 1
P=zeros(N+1,N+1,nlat);
for l=0:N
    P(l+1,1:l+1,:)=legendre(l,cos(theta),'norm');
end

damp=exp(-(0:N)'.*((0:N)'+1)*sigma^2/2); % heat kernel for each degree l

F=fft(data,[],2); % zonal wavenumbers m: column m+1 (and nlon-m+1 for -m)
G=zeros(size(F));
for m=0:min(N,floor((nlon-1)/2))
    Pm=reshape(P(m+1:N+1,m+1,:),N-m+1,nlat); % degrees l=m:N
    A=Pm'*(damp(m+1:N+1).*(Pm.*w')); % project, damp, synthesize
    G(:,m+1)=A*F(:,m+1);
    if m>0
        G(:,nlon-m+1)=A*F(:,nlon-m+1); % same filter for wavenumber -m
    end
end
outfield=real(ifft(G,[],2));

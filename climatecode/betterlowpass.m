function [tlow,ulowvect,uhivect]=betterlowpass(time,u,om,ns)

% [tlow,ulowvect,uhivect]=betterlowpass(time,u,om,ns)
% Lanczos lowpass filter applied to each column of a matrix, by calling
% lowpass.m (Emery and Thomson, Data Analysis Methods in Physical
% Oceanography, section 5.10).
%
% Inputs
% time = timestamps (vector of length size(u,1))
% u    = data, time down the rows, one series per column (evenly spaced)
% om   = 2*pi*dt*cutoff_frequency (= pi*cutoff/Nyquist), as in lowpass.m
% ns   = filter half-width in points (filter has 2*ns+1 weights)
%
% Outputs (size(u,1)-2*ns rows; ns points are lost at each end)
% tlow     = timestamps of filtered data
% ulowvect = lowpass filtered data, one column per column of u
% uhivect  = highpass filtered data (u - lowpass)
%
% NaN in u stays local: output within ns points of a NaN is NaN.
%
% Bug fixes 8 October 2026 (with Claude): output is no longer truncated
% at the first row where every column is NaN.  Previously a NaN gap
% present in all columns cut off all data after the gap; now the full
% size(u,1)-2*ns rows are returned.

[size5,size6]=size(u);

% lowpass drops ns points from each end of the series
nlow=max(size5-2*ns,0);
ulowvect=zeros(nlow,size6);
uhivect=zeros(nlow,size6);
for k=1:size6
    [tlow,ulowvect(:,k),uhivect(:,k)]=lowpass(time,u(:,k),om,ns);
end

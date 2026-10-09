function [EOFs,ECs,expv,varargout]=EOF(inmatrix,n,timedimension,weights,detrendflag,scaleflag,reconstructionsflag)

% Function EOF
% Calculates EOFs for a matrix via singular value decomposition
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 18 February 2020
% Bug fixes 8 October 2026 (with Claude): mean removal (detrendflag=1) now
%   subtracts the time mean at each spatial point (previously the spatial mean
%   at each time was removed, so EOFs were not of anomalies about the
%   climatology); eigenvalues lambda now use (number of times - 1), not (number
%   of spatial points - 1); expv is now the fraction of the TOTAL variance, so
%   it is correct when n>0 (previously it summed to 1 over the n returned
%   modes); scaleflag=1 now divides each EOF by the square root of its own
%   eigenvalue, and the reconstructions undo this only when scaleflag=1
%   (previously both used a wrongly shaped factor, the sum of sqrt(lambda),
%   so scaled EOFs and all reconstructions were wrong); latitude weighting (weights<0) works for
%   any latitude dimension and number of dimensions; weights, detrendflag,
%   scaleflag and reconstructionsflag are now optional; the weights size is
%   checked; economy-size SVD for n=0.  Grid points with weight 0 or with
%   any NaN are now excluded from the SVD and come out as NaN in the EOFs
%   and reconstructions (previously weight 0 gave NaN/Inf from dividing by
%   sqrt(weights), and any NaN broke the SVD); results are unchanged when
%   there are no zero weights or NaNs.
%
% [EOFs,ECs,expv,reconstructions]=EOF(inmatrix,n,timedimension,weights,detrendflag,scaleflag,reconstructionsflag)
%
% Inputs:
%   inmatrix = the input matrix (you do NOT need to transform it).  Any
%              number of dimensions; one of them is time, the others are
%              treated as space and flattened.  Any spatial point with a
%              NaN at any time is excluded from the analysis (as if its
%              weight were 0).
%   n        = The number of EOFs you want returned.  If n=0, then it does
%              them all (min(number of spatial points, number of times)).
%              If n>0, svds is used to compute only the leading n modes.
%   timedimension = the dimension of inmatrix that represents time
%   weights = (optional, default 1) A map of weights for your input
%             variables (rows).  For geospatial
%             data, this is often an area weighting map to deemphasize high
%             latitude grid points that may have high variance but do not contribute
%             much to patterns.  You can also do masking by setting certain
%             weights to 0: those points (and any with NaNs) are left out of
%             the SVD, which is done on the remaining points only, and come
%             out as NaN in the EOFs and reconstructions (output sizes are
%             unchanged).  A map must have one value per spatial point,
%             ordered like the non-time dimensions of inmatrix.  The SVD is
%             done on the field multiplied by sqrt(weights), so that the
%             covariance is weighted by weights.
%             You can set this to scalar values:
%             If weights==1, then the script will set this field to all ones.
%             If weights<0, then the script will calculate an area weighting
%                 map approximately based on the dimensions of your input matrix.
%                 It will assume the latitude dimension is the modulus of the
%                 argument.  (For example, weights=-2 means latitude is the
%                 second dimension.)  This assumes a global grid with equally
%                 spaced latitudes from pole to pole (weights are the areas of
%                 the latitude bands).
%             Methodology from Baldwin et al. (2009), doi:10.1175/2008JCLI2147.1
%   detrendflag:  (optional, default 1) If 0, the routine will not remove
%                 the mean for you (this assumes inmatrix is already
%                 anomalies).  If 1, the routine will subtract the time
%                 mean at every spatial point (it does not remove linear
%                 trends).  I strongly recommend using detrendflag=1
%                 unless you are sure of what you are doing.
%   scaleflag:  (optional, default 0) If 1, the script will divide each EOF
%               by the square root of its eigenvalue (lambda, below), so
%               the EOFs are dimensionless patterns and the amplitude is
%               carried by lambda.  Note that
%               this rescaling will be undone when computing the
%               reconstructions.  (This does not standardize the input; if
%               you want rows/variables with high variance not to dominate,
%               standardize inmatrix before calling.)
%   reconstructionsflag:  (optional, default 0) If 0, the output will just be the EOFs (spatial
%                 patterns) and ECs (time series for those patterns).  If
%                 this is a positive integer, it will output reconstructed
%                 timeseries incorporating the first n EOFs.  For example,
%                 if this flag is 3, it will return the cell array
%                 {EC1*EOF1,EC1*EOF1+EC2*EOF2,EC1*EOF1+EC2*EOF2+EC3*EOF3}.
%                 It is capped at the number of modes computed.
%
% Outputs:
%   This function uses the singular value decomposition O=U*S*V', where O
%   is your (weighted, mean-removed) input matrix, flattened to
%   (spatial points) x (times).  The EOFs (spatial patterns) are U*S, and the ECs
%   (the timeseries associated with each EOF) are V.  The singular values
%   (useful for determining the explained variance of EOFs) are given by S.
%
%   EOFs = the EOFs, (number of spatial points) x (number of modes), in the
%          units of inmatrix (unless scaleflag=1).  Spatial points are in
%          column-major order of the non-time dimensions, so use reshape to
%          turn a column back into a map.
%   ECs  = the ECs, (number of times) x (number of modes); each column has
%          unit length (norm 1, not unit variance)
%   expv = the fraction (0 to 1) of the total weighted variance explained
%          by each of the modes, (number of modes) x 1 (total over the
%          included points only, i.e. excluding masked/NaN points)
%   reconstructions (optional) = 1 x reconstructionsflag cell array; each
%          cell is the same size as inmatrix and includes the mean that was
%          removed (NaN at masked/NaN points).  Returned as {} if
%          reconstructionsflag=0.

% default values for optional inputs
if nargin<4 || isempty(weights)
    weights=1;
end
if nargin<5 || isempty(detrendflag)
    detrendflag=1;
end
if nargin<6 || isempty(scaleflag)
    scaleflag=0;
end
if nargin<7 || isempty(reconstructionsflag)
    reconstructionsflag=0;
end

% setting up the permutation matrices so that time is the last dimension when doing SVD
I=size(inmatrix);
J=length(I);
K=setdiff(1:J,timedimension);
permutematrix=[K timedimension];
if timedimension==1
    reversepermutematrix=[J 1:J-1];
else
    reversepermutematrix=[1:timedimension-1 J timedimension:J-1];
end

% making the input matrix SVD-friendly (flattening spatial dimensions)
M=permute(inmatrix,permutematrix);
Q=size(M);
R=length(Q);
N=double(reshape(M,[prod(Q(1:R-1)),Q(R)])); % SVD likes double precision inputs
NI=size(N);

% creating matrix of weights
if isscalar(weights)
    if weights==1
        weights=ones([Q(1:R-1) 1]);
    elseif weights<0 && abs(weights)~=timedimension && abs(weights)<=J
        latdim=find(permutematrix==abs(weights)); % latitude dimension after permuting
        latedges=linspace(-90,90,Q(latdim)+1); % assumes equally spaced, pole-to-pole latitudes
        colatedges=pi/2-latedges*pi/180;
        % area of each latitude band (difference of cos(colatitude)), broadcast along the other dimensions
        wsize=ones(1,max(R-1,2));
        wsize(latdim)=Q(latdim);
        weights=reshape(abs(diff(cos(colatedges))),wsize).*ones([Q(1:R-1) 1]);
    else
        error('EOF:weights','Scalar weights must be 1 or the negative of a non-time dimension of inmatrix.');
    end
end
if numel(weights)~=NI(1)
    error('EOF:weights','weights must have one value per spatial point (%d).',NI(1));
end

% excluding points with weight 0 or any NaN from the SVD (they are set to NaN in the outputs)
w=double(weights(:));
keep=(w~=0) & ~any(isnan(N),2);
if ~any(keep)
    error('EOF:mask','No spatial points left after excluding zero weights and NaNs.');
end
N=N(keep,:);

if detrendflag==1
    O=N-mean(N,2); % subtracts the time mean from every row (grid point)
    Pk=N-O; % saving the mean values so you can add them back in later
else
    O=N;
    Pk=zeros(size(N));
end
P=NaN(NI); % means on the full grid (NaN at excluded points)
P(keep,:)=Pk;

sqw=sqrt(w(keep)); % square root because O'*O should be weighted by weights

O=O.*sqw; % weighting the field on which we will perform SVD

if n==0 % do all svds
    [U,S,V]=svd(O,'econ');
else % do n svds
    [U,S,V]=svds(O,n);
end
m=size(S,2);

% creating outputs
EOFs=U*S;
ECs=V;
lambda=diag(S).^2/(NI(2)-1); % eigenvalues (variance of each EOF)
if scaleflag==1
    EOFs=EOFs./transpose(sqrt(lambda)); % dividing each EOF (column) by the square root of its eigenvalue
end
expv=diag(S).^2/sum(O(:).^2); % diag over trace (total variance of included points, including modes not returned)

% unweighting the EOFs and putting them back on the full grid (NaN at excluded points)
EOFk=EOFs./sqw;
EOFs=NaN(NI(1),m);
EOFs(keep,:)=EOFk;

% outputs the reconstructions (if requested)
if reconstructionsflag>0
    if reconstructionsflag>m
        reconstructionsflag=m;
    end
    outr=cell(1,reconstructionsflag);
    for k=1:reconstructionsflag
        if scaleflag==1
            recEOFs=EOFs(:,1:k).*transpose(sqrt(lambda(1:k))); % undoing the rescaling
        else
            recEOFs=EOFs(:,1:k);
        end
        outr{k}=permute(reshape(recEOFs*transpose(ECs(:,1:k))+P,Q),reversepermutematrix);
    end
    varargout{1}=outr;
elseif nargout>3
    varargout{1}={};
end

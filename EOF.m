function [EOFs,ECs,expv,varargout]=EOF(inmatrix,n,timedimension,weights,detrendflag,scaleflag,reconstructionsflag)

% Function EOF
% Calculates EOFs for a matrix via singular value decomposition
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 18 February 2020
%
% Inputs:
%   inmatrix = the input matrix (you do NOT need to transform it)
%   n        = The number of EOFs you want returned.  If n=0, then it does
%              them all (up to the number of times you have).
%   timedimension = the dimension of inmatrix that represents time
%   weights = A map of weights for your input variables (rows).  For geospatial
%             data, this is often an area weighting map to deemphasize high
%             latitude grid points that may have high variance but do not contribute
%             much to patterns.  You can also do masking by setting certain
%             weights to 0.  You can set this to scalar values:
%             If weights==1, then the script will set this field to all ones.
%             If weights<0, then the script will calculate an area weighting
%                 map approximately based on the dimensions of your input matrix.
%                 It will assume the latitude dimension is the modulus of the
%                 argument.  (For example, weights=-2 means latitude is the
%                 second dimension.)
%             Methodology from Baldwin et al. (2009), doi:10.1175/2008JCLI2147.1
%   detrendflag:  If 0, the routine will not detrend data for you (this
%                 assumes inmatrix has already been detrended).  If 1,
%                 the routine will subtract the mean from all spatial
%                 patterns.  I strongly recommend using detrendflag=1
%                 unless you are sure of what you are doing.
%   scaleflag:  If 1, the script will rescale all EOFs to have
%               integer variance.  This is recommended so that rows/variables
%               with high variance do not dominate the analysis.  Note that
%               this rescaling will be undone when computing the
%               reconstructions.
%   reconstructionsflag:  If 0, the output will just be the EOFs (spatial
%                 patterns) and ECs (time series for those patterns).  If
%                 this is a positive integer, it will output reconstructed
%                 timeseries incorporating the first n EOFs.  For example,
%                 if this flag is 3, it will return the cell array
%                 {EC1*EOF1,EC1*EOF1+EC2*EOF2,EC1*EOF1+EC2*EOF2+EC3*EOF3}.
%
% Outputs:
%   This function uses the singular value decomposition M=U*S*V, where M
%   is your input matrix.  The EOFs (spatial patterns) are U*S, and the ECs
%   (the timeseries associated with each EOF) are V.  The singular values
%   (useful for determining the explained variance of EOFs) are given by S.
%
%   EOFs = the EOFs
%   ECs  = the ECs
%   expv = the explained variance of each of the modes
%   reconstructions (optional)

% setting up the permutation matrices so that time is the last dimension when doing SVD
I=size(inmatrix);
J=length(I);
K=setdiff(1:J,timedimension);
permutematrix=[K timedimension];
if timedimension==1;
    reversepermutematrix=[J 1:J-1];
else
    reversepermutematrix=[1:timedimension-1 J timedimension:J-1];
end

% making the input matrix SVD-friendly (flattening spatial dimensions)
M=permute(inmatrix,permutematrix);
Q=size(M);
R=length(Q);
N=reshape(M,[prod(Q(1:R-1)),Q(R)]);
NI=size(N);

if detrendflag==1;
    O=detrend(N,'constant'); % subtracts the mean value from every column
    P=N-O; % saving the mean values so you can add them back in later
else
    O=N;
    P=zeros(size(N));
end

% creating matrix of weights
if length(weights)==1;
    if weights==1;
        weights=ones(Q(1:R-1));
    elseif weights<0;
        latedges=linspace(-90,90,Q(abs(weights))+1);
        colatedges=pi/2-latedges*pi/180;
        weights=(repmat(abs(diff(cos(colatedges))),[Q(1) 1]));
    end
end
weights=reshape(weights,[prod(Q(1:R-1)) 1]);
W12=diag(sqrt(weights));
W12inv=diag(1./sqrt(weights));

O=transpose(transpose(O)*W12); % weighting the field on which we will perform SVD
O=double(O); % SVD likes double precision inputs

if n==0; % do all svds
    [U,S,V]=svd(O);
    m=Q(end);
else % do n svds
    [U,S,V]=svds(O,n);
    m=n;
end

% creating outputs
EOFs=U*S;
ECs=V;
lambda=diag(S).^2/(NI(1)-1); % eigenvalues (variance of each EOF)
if scaleflag==1;
    EOFs=EOFs./(ones(size(EOFs))*sqrt(lambda)); % dividing by square root of eigenvalues
end
expv=lambda/(sum(lambda)); % diag over trace

% unweighting the EOFs
EOFs=W12inv*EOFs;

% outputs the reconstructions (if requested)
if reconstructionsflag>0;
    if reconstructionsflag>m;
        reconstructionsflag=m;
    end
    outr=cell(1,reconstructionsflag);
    for k=1:reconstructionsflag;
        outr{k}=permute(reshape(EOFs(:,1:k).*(ones(size(EOFs(:,1:k)))*sqrt(lambda(1:k,:)))*transpose(ECs(:,1:k))+P,Q),reversepermutematrix);
    end
    varargout{1}=outr;
end


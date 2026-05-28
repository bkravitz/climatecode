function outvals=tripolarconversion(inputfield,varname,inputlatvertices,inputlonvertices,latbnds,lonbnds)

% function outvals=tripolarconversion(filename,varname,latbnds,lonbnds)
% Written by Ben Kravitz (ben.kravitz@pnnl.gov or ben.kravitz.work@gmail.com) 9 May 2014
%
% This function is designed to convert tripolar ocean grids (or any grids characterized
% by irregular quadrilaterals) to lat-lon grids.
%
% Inputs:  inputfield       - the data you are converting
%          varname          - the variable name of the data
%          inputlatvertices - an array of the latitudes of the vertices of each quadrilateral
%          inputlonvertices - an array of the latitudes of the vertices of each quadrilateral
%          latbnds          - the boundaries (edges) of the latitudes that you want, e.g., [-90 -88.75-0.625:1.25:88.75+0.625 90]
%          lonbnds  - the boundaries (edges) of the longitudes that you want
%
% This script assumes that there are variables 'lat_vertices' and 'lon_vertices' in the input
%   file that describe the lat-lon coordinates of the vertices of the irregular quadrilaterals.
%   This should be the case for CMORized data, but it may not be true generically.  The script
%   should be easy to modify if you need to specify these vertices.

I=size(inputfield);
if numel(I)<3;
    I(3)=1;
end
% flattening the ocean array
tos2=zeros(I(1)*I(2),I(3));
for k=1:I(1);
    tos2((k-1)*I(2)+1:k*I(2),:)=inputfield(k,:,:);
end

% flattening the arrays of the grid "box" vertices
latv2=zeros(4,I(1)*I(2));
lonv2=zeros(4,I(1)*I(2));
for k=1:I(1);
    latv2(:,(k-1)*I(2)+1:k*I(2))=inputlatvertices(:,k,:);
    lonv2(:,(k-1)*I(2)+1:k*I(2))=inputlonvertices(:,k,:);
end

% fixing the longitude list in case a quadrilateral spans the 0 longitude meridian
for k=1:length(lonv2);
    blah=lonv2(:,k);
    a=[abs(blah(1)-blah(2)) abs(blah(1)-blah(3)) abs(blah(1)-blah(4)) abs(blah(2)-blah(3)) abs(blah(2)-blah(4)) abs(blah(3)-blah(4))];
    Q=find(a>300);
    if isempty(Q)==0;
        J=find(blah<180);
        blah(J)=blah(J)+360;
        lonv2(:,k)=blah;
    end
end

% the latitude and longitude bounds that you actually want
II=size(latbnds);
if II(2)>II(1);latbnds=transpose(latbnds);end
II=size(lonbnds);
if II(2)>II(1);lonbnds=transpose(lonbnds);end
lat2=transpose(repmat(latbnds,[1 length(lonbnds)]));
lon2=repmat(lonbnds,[1 length(latbnds)]);
J=length(latbnds);
K=length(lonbnds);

% formulating a Delaunay Triangulation of the domain
% Essentially, what we are doing is:
%   1) Figure out all of the vertices of the tripolar grid
%   2) Figure out all of the vertices of the lat-lon grid
%   3) Collect all of these vertices into one huge mass of points
%   4) Form a triangulation with all of these points such that no triangle overlaps.  This is called a Delaunay Triangulation.
%      The advantage is that now each triangle fits exactly within one irregular quadrilateral from the polar grid and exactly
%      one lat-lon grid box.
%   5) Then we assign a value (e.g., sea surface temperature) to each triangle based on the value of the irregular quadrilateral.
%      Inside each lat-lon grid box, we then perform an area-weighted average of all of the triangles contained in that box.
%   6) Then we win.

% Forming the Delaunay Triangulation using a built in matlab function
% The triangulation will output a structure with two entries:
% DTX.X is a list of the coordinates of the vertices (in lat-lon)
% DTX.Triangulation gives a list of the vertices associated with each simplex (triangle).  The simplices are numbered by matlab.
DTX=DelaunayTri(cat(1,double(reshape(inputlonvertices,[4*I(1)*I(2) 1])),reshape(lon2,[K*J 1])),cat(1,double(reshape(inputlatvertices,[4*I(1)*I(2) 1])),reshape(lat2,[K*J 1])));
SV=DTX.Triangulation;

% returns the centers of the simplices (triangles) created by the triangulation (in lat-lon)
Xcenters=incenters(DTX);

% This finds the area of each simplex using a formula for an arbitrary triangle
% This is done using the determinant function (det) and can take a bit of time.
areass=zeros(length(Xcenters),1);
for k=1:length(Xcenters);
    areass(k)=0.5*abs(det([DTX.X(SV(k,1),1) DTX.X(SV(k,2),1) DTX.X(SV(k,3),1) ; DTX.X(SV(k,1),2) DTX.X(SV(k,2),2) DTX.X(SV(k,3),2) ; 1 1 1]*pi/180));
end

% if a simplex is in a grid "box" (determined using inpolygon), give the simplex that value (vallist)
% inpolygon can take a bit of time
% if you are not careful, this step can use a LOT of memory
vallist=zeros(length(Xcenters),I(3));
valcount=zeros(length(Xcenters),I(3));
for j=1:length(lonv2);
    INC=inpolygon(Xcenters(:,1),Xcenters(:,2),[lonv2(1,j) ; lonv2(2,j) ; lonv2(3,j) ; lonv2(4,j) ;lonv2(1,j)],[latv2(1,j) ; latv2(2,j) ; latv2(3,j) ; latv2(4,j) ;latv2(1,j)]);
    Q=find(INC==1);
    JJ=1-isnan(tos2(j,:));
    blah=squeeze(nansum(cat(3,vallist(Q,:),repmat(tos2(j,:),[length(Q) 1])),3));
    vallist(Q,:)=blah;
    blah=valcount(Q,:)+repmat(JJ,[length(Q) 1]).*ones(length(Q),I(3));
    valcount(Q,:)=blah;
end

% determining the average value in each simplex
% this is sort of a formality because the value in each simplex should be defined by the "box" containing it, but this is a clean way of taking care of missing values
vallist=vallist./valcount;

outvals=zeros(K-1,J-1,I(3));
outcount=zeros(K-1,J-1);
for i=1:K-1;
    for j=1:J-1;
        IN=inpolygon(Xcenters(:,1),Xcenters(:,2),[lonbnds(i) lonbnds(i+1) lonbnds(i+1) lonbnds(i) lonbnds(i)],[latbnds(j) latbnds(j) latbnds(j+1) latbnds(j+1) latbnds(j)]);
        Q=find(IN==1);
        QQ=find(sum(isnan(vallist(Q,:))==1,2)>0);
        IN(Q(QQ))=0;
        Q=find(IN==1);
        outcount(i,j)=sum(areass.*IN);
        outvals(i,j,:)=outvals(i,j,:)+reshape(nansum(vallist(Q,:).*repmat(areass(Q).*IN(Q),[1 I(3)]),1),[1 1 I(3)]);
    end
end
outcount=repmat(outcount,[1 1 I(3)]);
outvals=outvals./outcount;

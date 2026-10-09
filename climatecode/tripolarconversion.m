function outvals=tripolarconversion(inputfield,varname,inputlatvertices,inputlonvertices,latbnds,lonbnds) %#ok<INUSL>

% TRIPOLARCONVERSION  Area-weighted regridding of quadrilateral (e.g. tripolar ocean) grids to lat-lon
%
% outvals=tripolarconversion(inputfield,varname,inputlatvertices,inputlonvertices,latbnds,lonbnds)
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com) 9 May 2014
%
% This function is designed to convert tripolar ocean grids (or any grids characterized
% by irregular quadrilaterals) to lat-lon grids.  Each lat-lon box gets the
% area-weighted mean of the quadrilaterals that overlap it (first-order
% conservative remapping).  Missing values (NaN, e.g. land) are left out of
% the mean; boxes with no valid data are NaN.
%
% Inputs:  inputfield       - the data you are converting, I1 x I2 x nt (nt optional, e.g. time)
%          varname          - the variable name of the data (not used; kept for compatibility)
%          inputlatvertices - the latitudes (degrees) of the 4 vertices of each quadrilateral,
%                             4 x I1 x I2 (I1 x I2 x 4 is also accepted), vertices in
%                             order around the cell (as in CMORized lat_vertices/vertices_latitude)
%          inputlonvertices - the longitudes (degrees) of the vertices, same size; any
%                             convention (0-360 or -180-180)
%          latbnds          - the boundaries (edges) of the latitudes that you want, e.g., [-90 -88.75-0.625:1.25:88.75+0.625 90]
%          lonbnds          - the boundaries (edges) of the longitudes that you want, e.g. 0:2.5:360
%                             (any convention; the input longitudes are wrapped to match)
%
% Output:  outvals          - (numel(lonbnds)-1) x (numel(latbnds)-1) x nt, i.e. lon x lat x nt
%
% Method:
%   1) Collect the vertices of the quadrilaterals and of the lat-lon grid.
%   2) Form a constrained Delaunay triangulation of all of these points, with
%      the edges of both grids as constraints (intersections of edges are added
%      as points).  Each triangle then lies in exactly one quadrilateral and
%      at most one lat-lon box.
%   3) Each triangle takes the value of the quadrilateral containing it.  Its
%      area is the planar area in lon-lat times cos(latitude) of its centroid
%      (i.e. the area on the sphere, to second order in the triangle size).
%   4) Inside each lat-lon box, take the area-weighted average of all of the
%      triangles contained in that box.
%   5) Then we win.
% Assumptions: cell edges are straight lines in lon-lat space, no
% quadrilateral contains a pole, and cells are not wider than 180 degrees
% of longitude.  A cell crossing the longitude seam of the output grid is
% split between the two ends of the grid.
%
% Bug fixes 8 October 2026 (with Claude):
%   - The Delaunay triangulation was unconstrained, so triangles could
%     straddle cell and box edges and be assigned to the wrong cell/box;
%     now constrained by both grids' edges.
%   - Triangle areas now include cos(latitude) (previously planar lon-lat
%     area, overweighting high latitudes within a box).
%   - Cells crossing the 0/360 meridian lost their part east of 0, and
%     input/output longitude conventions (0-360 vs -180-180) were not
%     reconciled; both fixed.
%   - A triangle missing data at any time was dropped at all times; now
%     the missing-data mask is applied separately for each time.
%   - Replaced DelaunayTri/incenters (deprecated) and nansum (Statistics
%     Toolbox) with base MATLAB; much faster (no inpolygon loops).
%   Results change near cell edges, near the seam, and slightly at high
%   latitudes.  Constant fields are reproduced exactly.

% ---- checking inputs and flattening the arrays
I=size(inputfield,[1 2]);
nt=numel(inputfield)/prod(I);
ncell=prod(I);
if ~isequal(size(inputlatvertices),size(inputlonvertices))
    error('tripolarconversion: lat and lon vertex arrays must be the same size.');
end
if size(inputlatvertices,1)~=4 && size(inputlatvertices,3)==4 % I1 x I2 x 4
    inputlatvertices=permute(inputlatvertices,[3 1 2]);
    inputlonvertices=permute(inputlonvertices,[3 1 2]);
end
if numel(inputlatvertices)~=4*ncell || size(inputlatvertices,1)~=4
    error('tripolarconversion: vertex arrays must be 4 x I1 x I2, matching inputfield.');
end
latbnds=double(latbnds(:));
lonbnds=double(lonbnds(:));
if numel(latbnds)<2 || numel(lonbnds)<2 || any(diff(latbnds)<=0) || any(diff(lonbnds)<=0)
    error('tripolarconversion: latbnds and lonbnds must be increasing, with at least 2 values.');
end
if lonbnds(end)-lonbnds(1)>360+1e-9
    error('tripolarconversion: lonbnds must span at most 360 degrees.');
end
J=numel(latbnds);
K=numel(lonbnds);

vals=reshape(double(inputfield),ncell,nt);          % cell values, ncell x nt
latv2=reshape(double(inputlatvertices),4,ncell);    % vertices, 4 x ncell
lonv2=reshape(double(inputlonvertices),4,ncell);

% ---- longitudes: wrap to [lonbnds(1), lonbnds(1)+360) to match the output
% grid, then unwrap cells that cross that seam and add a copy of each of them
% shifted by -360, so both ends of the output grid are covered
lon0=lonbnds(1);
lonv2=lon0+mod(lonv2-lon0,360);
cross=(max(lonv2,[],1)-min(lonv2,[],1))>180;
lonv2(:,cross)=lonv2(:,cross)+360*(lonv2(:,cross)<lon0+180);
cellid=[1:ncell find(cross)];                       % source cell of each polygon
lonv2=[lonv2 lonv2(:,cross)-360];
latv2=[latv2 latv2(:,cross)];
npoly=numel(cellid);

% ---- points and constraint edges of both grids
[LO,LA]=ndgrid(lonbnds,latbnds);                    % lat-lon grid nodes, K x J
nodeid=reshape(1:K*J,K,J);
gridedges=[reshape(nodeid(1:K-1,:),[],1) reshape(nodeid(2:K,:),[],1) ; ...  % along latitude lines
    reshape(nodeid(:,1:J-1),[],1) reshape(nodeid(:,2:J),[],1)];             % along meridians
pts=[lonv2(:) latv2(:) ; LO(:) LA(:)];
quadid=reshape(1:4*npoly,4,npoly);
quadedges=[reshape(quadid,[],1) reshape(quadid([2 3 4 1],:),[],1)];
edges=[quadedges ; gridedges+4*npoly];

% remove duplicate points (shared vertices) and degenerate/duplicate edges
[P,~,ic]=uniquetol(pts,1e-10,'ByRows',true,'DataScale',1);
edges=ic(edges);
edges=edges(edges(:,1)~=edges(:,2),:);
edges=unique(sort(edges,2),'rows');

% ---- constrained Delaunay triangulation
% (splitting intersecting edges is expected, so those warnings are turned off)
ws=[warning('off','MATLAB:delaunayTriangulation:ConsConsSplitWarnId') ...
    warning('off','MATLAB:delaunayTriangulation:ConsSplitPtWarnId') ...
    warning('off','MATLAB:delaunayTriangulation:DupPtsConsUpdatedWarnId')];
DT=delaunayTriangulation(P,edges);
warning(ws);
T=DT.ConnectivityList;
X=DT.Points;
xc=(X(T(:,1),1)+X(T(:,2),1)+X(T(:,3),1))/3;         % triangle centroids (lon, lat)
yc=(X(T(:,1),2)+X(T(:,2),2)+X(T(:,3),2))/3;
% area of each triangle on the unit sphere: planar area (radians^2) times cos(lat)
areass=0.5*abs((X(T(:,2),1)-X(T(:,1),1)).*(X(T(:,3),2)-X(T(:,1),2)) - ...
    (X(T(:,3),1)-X(T(:,1),1)).*(X(T(:,2),2)-X(T(:,1),2)))*(pi/180)^2.*cosd(yc);

% ---- which lat-lon box contains each triangle
ib=discretize(xc,lonbnds);
jb=discretize(yc,latbnds);

% ---- which quadrilateral contains each triangle: split each quadrilateral
% into two triangles along a diagonal that lies inside it, and locate the
% centroids in that triangulation
qx=lonv2;
qy=latv2;
% vertices 2 and 4 on opposite sides of diagonal 1-3 -> use 1-3, else 2-4
s2=(qx(3,:)-qx(1,:)).*(qy(2,:)-qy(1,:))-(qy(3,:)-qy(1,:)).*(qx(2,:)-qx(1,:));
s4=(qx(3,:)-qx(1,:)).*(qy(4,:)-qy(1,:))-(qy(3,:)-qy(1,:)).*(qx(4,:)-qx(1,:));
use13=(s2.*s4<=0);
tri1=[quadid(1,:) ; quadid(2,:) ; quadid(3,:)];
tri2=[quadid(1,:) ; quadid(3,:) ; quadid(4,:)];
tri1(:,~use13)=[quadid(2,~use13) ; quadid(3,~use13) ; quadid(4,~use13)];
tri2(:,~use13)=[quadid(2,~use13) ; quadid(4,~use13) ; quadid(1,~use13)];
TRq=triangulation([tri1' ; tri2'],[lonv2(:) latv2(:)]);
tq=pointLocation(TRq,[xc yc]);
src=NaN(size(tq));
ok=~isnan(tq);
src(ok)=cellid(mod(tq(ok)-1,npoly)+1);              % triangle -> polygon -> cell

% ---- area-weighted averages in each lat-lon box
use=~isnan(ib) & ~isnan(jb) & ~isnan(src) & areass>0;
box=sub2ind([K-1 J-1],ib(use),jb(use));
A=sparse(box,src(use),areass(use),(K-1)*(J-1),ncell); % overlap area of each box with each cell
valid=~isnan(vals);
vals(~valid)=0;
outvals=full(A*vals)./full(A*double(valid));          % 0/0 = NaN where no valid data
outvals=reshape(outvals,[K-1 J-1 nt]);

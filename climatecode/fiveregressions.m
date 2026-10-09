function [x,y,m,mvar,b,bvar,r]=fiveregressions(x,y)

% [x,y,m,mvar,b,bvar,r]=fiveregressions(x,y)
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 17 February 2020
% Bug fixes 8 October 2026 (with Claude): row-vector inputs now work
%   (previously centering produced n x n matrices and wrong results); the
%   orthogonal-regression intercept variance (bvar(4)) no longer uses
%   abs(m1), which gave the wrong value when the slope is negative; x and y
%   of different sizes raise an error.
%
% This program calculates regressions five different ways, depending upon
% the independence/dependence of the input variables x and y.
% The function will output for each of the five regression methods
% slope (m), y-intercepts (b), and the variances of each (mvar and bvar;
% take sqrt for one standard deviation error bars).  The function will
% also output r, the correlation coefficient.
%
% Method and variance formulas (asymptotic, valid for large n) follow
% Isobe et al. (1990), Linear regression in astronomy I, ApJ 364, 104.
%
% Inputs:  x, y = vectors (or arrays) with the same number of elements;
%                 pairs where either is NaN are removed
% Outputs: x, y = the NaN-free data actually used, as column vectors
%          m, mvar, b, bvar = 1 x 5 rows, one entry per method (order below);
%                 all slopes are expressed as dy/dx
%          r    = Pearson correlation coefficient (scalar)
%
% The five regression methods are (in order):
% (1) Ordinary least squares regression with x as the independent variable
%     and y as the dependent variable
% (2) Ordinary least squares regression with y as the independent variable
%     and x as the dependent variable (slope reported as dy/dx)
% (3) Ordinary least squares bisector method, which bisects methods 1 and 2
% (4) Orthogonal regression (total least squares regression), an
%     errors-in-variables method
% (5) Reduced major axis (RMA) regression, an errors-in-variables method

if numel(x)~=numel(y)
    error('fiveregressions:sizeMismatch','x and y must have the same number of elements');
end

% Removing all NaN values from x and y (pairwise)
x=x(:);
y=y(:);
I=~isnan(x) & ~isnan(y);
x=x(I);
y=y(I);
n=length(x);

% centering x and y and computing some commonly used values
xbar=sum(x)/n;
ybar=sum(y)/n;
x0=x-xbar;
y0=y-ybar;
sxx=sum(x0.^2);
syy=sum(y0.^2);
sxy=sum(x0.*y0);

% formulas for slopes for each regression method
m1=sxy/sxx;
m2=syy/sxy;
m3=(m1*m2-1+sqrt((1+m1^2)*(1+m2^2)))/(m1+m2);
%m4=0.5*(m2-1/m1)+sign(sxy)*sqrt(4+(m2-1/m1)^2);
m4=(syy-sxx+sqrt((syy-sxx)^2+4*sxy^2))/(2*sxy);
m5=sign(sxy)*sqrt(m1*m2);

% slope variances
m1var=(1/sxx.^2).*sum(x0.^2.*(y-m1*x-ybar+m1*xbar).^2);
m2var=(1/sxy.^2).*sum(y0.^2.*(y-m2*x-ybar+m2*xbar).^2);
covm1m2=sum(x0.*y0.*(y0-m1*x0).*(y0-m2*x0))/(m1*sxx^2);
m3var=(m3^2/((m1+m2)^2*(1+m1^2)*(1+m2^2)))*((1+m2^2)^2*m1var+2*(1+m1^2)*(1+m2^2)*covm1m2+(1+m1^2)^2*m2var);
m4var=(m4^2/(4*m1^2+(m1*m2-1)^2))*(m1var/m1^2+2*covm1m2+m2var*m1^2);
m5var=0.25*(m2*m1var/m1+2*covm1m2+m1*m2var/m2);

% y-intercepts
b1=ybar-m1*xbar;
b2=ybar-m2*xbar;
b3=ybar-m3*xbar;
b4=ybar-m4*xbar;
b5=ybar-m5*xbar;

% computing some useful quantities for y-intercept variances
% (gamma1, gamma2 = partial derivatives of each slope with respect to m1, m2)
g1=m3/((m1+m2)*sqrt((1+m1^2)*(1+m2^2)));
g2=m4/sqrt(4*m1^2+(m1*m2-1)^2);
gamma1=[1 0 g1*(1+m2^2) g2/m1 0.5*sqrt(m2/m1)];
gamma2=[0 1 g1*(1+m1^2) g2*m1 0.5*sqrt(m1/m2)];

% y-intercept variances
b1var=(1/n^2)*sum((y0-m1*x0-n*xbar*((gamma1(1)/sxx)*x0.*(y0-m1*x0)+(gamma2(1)/sxy)*y0.*(y0-m2*x0))).^2);
b2var=(1/n^2)*sum((y0-m2*x0-n*xbar*((gamma1(2)/sxx)*x0.*(y0-m1*x0)+(gamma2(2)/sxy)*y0.*(y0-m2*x0))).^2);
b3var=(1/n^2)*sum((y0-m3*x0-n*xbar*((gamma1(3)/sxx)*x0.*(y0-m1*x0)+(gamma2(3)/sxy)*y0.*(y0-m2*x0))).^2);
b4var=(1/n^2)*sum((y0-m4*x0-n*xbar*((gamma1(4)/sxx)*x0.*(y0-m1*x0)+(gamma2(4)/sxy)*y0.*(y0-m2*x0))).^2);
b5var=(1/n^2)*sum((y0-m5*x0-n*xbar*((gamma1(5)/sxx)*x0.*(y0-m1*x0)+(gamma2(5)/sxy)*y0.*(y0-m2*x0))).^2);

% outputs
m=[m1 m2 m3 m4 m5];
mvar=[m1var m2var m3var m4var m5var];
b=[b1 b2 b3 b4 b5];
bvar=[b1var b2var b3var b4var b5var];
r=sxy/sqrt(syy*sxx);

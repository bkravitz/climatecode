function [x,y,m,mvar,b,bvar,r]=fiveregressions(x,y)

% [x,y,m,mvar,b,bvar,r]=fiveregressions(x,y)
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 17 February 2020
%
% This program calculates regressions five different ways, depending upon
% the independence/dependence of the input variables x and y.
% The function will output for each of the five regression methods
% slope (m), y-intercepts (b), and error bars in each (mvar and bvar),
% represented as one standard deviation.  The function will also output
% r, the correlation coefficient.
%
% The five regression methods are (in order):
% (1) Ordinary least squares regression with x as the independent variable
%     and y as the dependent variable
% (2) Ordinary least squares regression with y as the independent variable
%     and x as the dependent variable
% (3) Ordinary least squares bisector method, which bisects methods 1 and 2
% (4) Orthogonal regression (total least squares regression), an
%     errors-in-variables method
% (5) Reduced major axis (RMA) regression, an errors-in-variables method

% Removing all NaN values from x and y
I=find(isnan(x)==0);
x=x(I);
y=y(I);
J=find(isnan(y)==0);
x=x(J);
y=y(J);
n=length(x);

% centering x and y and computing some commonly used values
xbar=sum(x)/n;
ybar=sum(y)/n;
x0=x-repmat(xbar,[n 1]);
y0=y-repmat(ybar,[n 1]);
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
m1var=(1/sxx.^2).*sum(x0.^2.*(rad-m1*tas-ybar+m1*xbar).^2);
m2var=(1/sxy.^2).*sum(y0.^2.*(rad-m2*tas-ybar+m2*xbar).^2);
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
g1=m3/((m1+m2)*sqrt((1+m1^2)*(1+m2^2)));
g2=m4/sqrt(4*m1^2+(m1*m2-1)^2);
gamma1=[1 0 g1*(1+m2^2) g2/abs(m1) 0.5*sqrt(m2/m1)];
gamma2=[0 1 g1*(1+m1^2) g2*abs(m1) 0.5*sqrt(m1/m2)];

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
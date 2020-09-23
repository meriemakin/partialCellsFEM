%derivative of analytical solution at one single point
function [dudx,dudy] = derivativeAnalyticalsolution(A,x,y,xC,yC,bc2,radius2);

r = sqrt((xC- x)^2.0 + (yC - y)^2.0);
drdx = (x - xC)*(r^2)^(-0.5);
drdy = (y - yC)*(r^2)^(-0.5);
c1 = -(radius2*bc2)/A;
dudx = c1 * (1/r) * drdx;
dudy = c1 * (1/r) * drdy;
%dudx = ((-radius2*bc2)/A)*((1/r)*drdx);
%dudy = ((-radius2*bc2)/A)*((1/r)*drdy);

%derivative of analytical solution at one single point
function [dudx,dudy] = derivativeAnalyticalsolution(x,y,xC,yC,bc1,bc2,radius1, radius2);

r = sqrt((xC- x)^2.0 + (yC - y)^2.0);
drdx = (x - xC)*(r^2)^(-0.5);
drdy = (y - yC)*(r^2)^(-0.5);
c1  = (bc1 - bc2)/(log(radius1) - log(radius2));
dudx = c1*((1/r)*drdx);
dudy = c1*((1/r)*drdy);

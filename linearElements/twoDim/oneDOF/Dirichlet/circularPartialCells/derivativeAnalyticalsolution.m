%derivative of analytical solution at one single point
function [dudx,dudy] = derivativeAnalyticalsolution(f,A,x,y,xC,yC);

dudx = -(f/2*A)*(x-xC);
dudy = -(f/2*A)*(y-yC);


%derivative of analytical solution at one single point
function [dudx,dudy] = derivativeAnalyticalsolution(f,A,x,y,xC,yC);


dudx = -(f/A)*(x-xC);
dudy = -(f/A)*(y-yC);

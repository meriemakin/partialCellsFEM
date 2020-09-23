function [u] = analyticalsolutionSinglePoint(bc1,bc2,x,y,xC,yC,radius1, radius2);

r  = sqrt((xC- x)^2.0 + (yC - y)^2.0);
c1 = (bc1 - bc2)/(log(radius1) - log(radius2)); 
c2 = bc1 - c1 * log(radius1);
u  = c1*log(r)+c2;  

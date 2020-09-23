function [u] = analyticalsolutionSinglePoint(f,A,bc1,bc2,x,y,xC,yC,radius,H);

r = sqrt((x-xC)^2 + (y-yC)^2);

u = bc1 - (f/(4*A))*(r^2.0 - radius^2.0) + (bc2/H) + ((f*radius)/ (2*H));
%u = ((A*radius)/(2*H)) + (K/H)*bc2 + bc1 + (A/(4*K))*(radius^2 - r^2);

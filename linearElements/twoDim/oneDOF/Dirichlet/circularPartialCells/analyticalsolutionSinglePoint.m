function [u] = analyticalsolutionSinglePoint(f,A,bc,x,y,xC,yC,radius);

u = bc + (f/4*A)*radius^2 - (f/4*A)*((x-xC)^2 + (y-yC)^2);

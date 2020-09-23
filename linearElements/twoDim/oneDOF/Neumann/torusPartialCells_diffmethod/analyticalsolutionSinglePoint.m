function [u] = analyticalsolutionSinglePoint(A,bc1,bc2,x,y,xC,yC,radius1, radius2,rhs);

r = sqrt((xC- x)^2.0 + (yC - y)^2.0);
c1 = -((radius2 * bc2)/A);
c2 = bc1 - c1 * log(radius1) + (rhs/(4*A*radius1^2));
u = c1 * log(r) + c2;
%u = bc1 + (((-radius2*bc2)/A)*(log(r)-log(radius1)));

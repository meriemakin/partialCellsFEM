%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius1 is the radius of the circular physical domain
function [y] = analyticalsolution(A,x_n,y_n,xC,yC,radius1,radius2,bc1,bc2,rhs);

n = length(x_n);
y = [];

c1 = -((radius2*bc2)/(A)) ;
c2 = bc1 - c1* log(radius1) ;
for i=1:n
    r = sqrt((xC- x_n(i))^2.0 + (yC - y_n(i))^2.0);
    if((r <= radius1) && (r >= radius2))
       sol = c1 * log(r) + c2;
       %sol = bc1 + (((-radius2*bc2)/A)*(log(r)-log(radius1)));
       y = [y;x_n(i), y_n(i),sol];
    end
end

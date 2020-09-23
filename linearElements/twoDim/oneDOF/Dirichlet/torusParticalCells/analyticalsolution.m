%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius1 is the radius of the circular physical domain
function [y] = analyticalsolution(x_n,y_n,xC,yC,radius1,radius2,bc1,bc2);

n = length(x_n);
y = [];

pos = 1;

c1  = (bc1 - bc2)/(log(radius1) - log(radius2));
c2  = bc1 - (c1 * log(radius1));

for i=1:n
    r = sqrt((xC- x_n(i))^2.0 + (yC - y_n(i))^2.0);
    if((r <= radius1) && (r >= radius2))
       sol = c1 * log(r) + c2;
       y = [y;x_n(i), y_n(i),sol];
    end
end

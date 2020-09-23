%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [y] = analyticalsolution(f,A,x_n,y_n,xC,yC,radius,bc1,bc2);

n = length(x_n);
y = [];

pos = 1;

for i=1:n
    r = sqrt((xC- x_n(i))^2.0 + (yC - y_n(i))^2.0);
    if(r <= radius)
       %sol = ((A*radius)/(2*H)) + (K/H)*bc2 + bc1 + (A/(4*K))*(radius^2 - r^2);
       sol = bc1 - (f/(4*A))*r^2.0;
       y = [y;x_n(i), y_n(i),sol];
       %pos = pos + 1;
    end
end

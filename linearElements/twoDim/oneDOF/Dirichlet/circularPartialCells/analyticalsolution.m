%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [y] = analyticalsolution(f,A,x_n,y_n,xC,yC,radius,bc);

n = length(x_n);
y = [];

pos = 1;

for i=1:n
    if(sqrt((xC- x_n(i))^2.0 + (yC - y_n(i))^2.0) <= radius)
       sol = bc + ((f/(4*A)))*radius^2 - ((f/(4*A))*((x_n(i)-xC)^2 + (y_n(i)-yC)^2));
       y = [y;x_n(i), y_n(i),sol];
       %pos = pos + 1;
    end
end

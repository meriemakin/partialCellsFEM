%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [y] = analyticalsolution(f,A,x_n,y_n,xC,yC,radius,bc);

n = length(x_n);
y = [];


for i=1:n
       sol = bc + ((f/(4*A)))*radius^2 - ((f/(4*A))*((x_n(i)-xC)^2 + (y_n(i)-yC)^2));
       y = [y;sol];
    end
end

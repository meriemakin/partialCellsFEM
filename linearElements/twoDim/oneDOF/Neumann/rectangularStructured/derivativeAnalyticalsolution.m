%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [dudx,dudy] = derivativeAnalyticalsolution(x,y,a,b,f,A,n);


dudx = 0;
dudy = 0;


for m=1:n
    for k=1:n
         dudx = dudx + (((1-cos(k*pi))*(1-cos(m*pi))*(cos((k*pi*x)/a))*(sin((m*pi*y)/b))*((k*pi)/a))/(m*k*((m^2/b^2)+(k^2/a^2))));
         dudy = dudy + (((1-cos(k*pi))*(1-cos(m*pi))*(sin((k*pi*x)/a))*(cos((m*pi*y)/b))*((m*pi)/b))/(m*k*((m^2/b^2)+(k^2/a^2))));
    end
end

dudx = ((4*f)/(A*pi^4)) * dudx;
dudy = ((4*f)/(A*pi^4)) * dudy;

            


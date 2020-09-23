%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [anasol] = analyticalsolutionSinglePoint(x,y,a,b,f,A,n);

anasol = 0;

for m=1:n
    for k=1:n
         anasol = anasol +  ((1-cos(k*pi))*(1-cos(m*pi))*(sin((k*pi*x)/a))*(sin((m*pi*y)/b)))/(m*k*((m^2/b^2)+(k^2/a^2)));
    end
end

anasol = ((4*f)/(A*pi^4)) * anasol;

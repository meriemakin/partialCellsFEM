%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [anasol] = analyticalsolution_copy(x_n,y_n,bc);

n = 300;
a = max(x_n);
b = max(y_n);

anasol = zeros(length(x_n),1);

for i =1:n

    l_i = b*sinh((i*pi*a)/b);
    m_i = a*sinh((i*pi*b)/a);
    A_i = (2*b/(i*pi*l_i))*bc*(1-cos(i*pi));
    B_i = (2*a/(i*pi*m_i))*bc*(1-cos(i*pi));

    for j=1:length(x_n) 
        x = x_n(j);
        y = y_n(j);
        anasol(j,1) = anasol(j,1) + A_i*sinh(((i*pi)/b)*(a-x))*sin(((i*pi)/b)*y) + A_i*sinh(((i*pi)/b)*(x))*sin(((i*pi)/b)*y) + B_i*sinh(((i*pi)/a)*(b-y))*sin(((i*pi)/a)*x) + B_i*sinh(((i*pi)/a)*(y))*sin(((i*pi)/a)*x);
    end

end

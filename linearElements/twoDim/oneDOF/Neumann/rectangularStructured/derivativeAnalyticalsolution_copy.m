%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC % 
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [dudx,dudy]] = derivativeAnalyticalsolution_copy(x,y,a,b,bc);

n = 200;

dudx = 0;
dudy = 0;

for i =1:n

    l_i = b*sinh((i*pi*a)/b);
    m_i = a*sinh((i*pi*b)/a);
    A_i = (2*b/(i*pi*l_i))*bc*(1-cos(i*pi));
    B_i = (2*a/(i*pi*m_i))*bc*(1-cos(i*pi));

    dudx = dudx + A_i*cosh(((i*pi)/b)*(a-x))*sin(((i*pi)/b)*y)*((-i*pi)/b)+ A_i*cosh(((i*pi)/b)*(x))*sin(((i*pi)/b)*y)*((i*pi)/b) + B_i*sinh(((i*pi)/a)*(b-y))*cos(((i*pi)/a)*x)*((i*pi)/a)) + B_i*sinh(((i*pi)/a)*(y))*cos(((i*pi)/a)*x)*((i*pi)/a);

    dudy = dudy + A_i*sinh(((i*pi)/b)*(a-x))*cos(((i*pi)/b)*y)*((i*pi)/b)+ A_i*sinh(((i*pi)/b)*(x))*cos(((i*pi)/b)*y)*((i*pi)/b) + B_i*cosh(((i*pi)/a)*(b-y))*sin(((i*pi)/a)*x)*((-i*pi)/a)) + B_i*cosh(((i*pi)/a)*(y))*sin(((i*pi)/a)*x)*((i*pi)/a);


end



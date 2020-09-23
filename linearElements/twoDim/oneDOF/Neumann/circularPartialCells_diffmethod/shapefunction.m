%---------------------------------------------------------------------------%                          
%                           Definition of the shape functions               %
%                           in the master element.                          %
%---------------------------------------------------------------------------%

%phi are the shape functions 
%xi and eta are the isoparametric coordinates
%nn is the number of nodes per element
function phi = shapefunction(xi,eta);

phi(1) = (1/4)*(1-xi)*(1-eta);
%phi(2) = (1/4)*(1+xi)*(1-eta);
phi(2) = (1/4)*(1-xi)*(1+eta);
phi(3) = (1/4)*(1+xi)*(1+eta);
%phi(4) = (1/4)*(1-xi)*(1+eta);
phi(4) = (1/4)*(1+xi)*(1-eta);

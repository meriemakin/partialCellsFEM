%---------------------------------------------------------------------------%                          
%                           Definition of the shape function derivatives    %
%                           in the master element.                          %
%---------------------------------------------------------------------------%

%Dphi_j is the derivative of the shape function  w.r.t. 
%the isoparametric coordinate.
%xi and eta are the isoparametric coordinates
function Dphi = shapefunctionFirstDerivative(xi,eta);

Dphi(1,1) = (-1/4)*(1-eta);
Dphi(1,2) = (-1/4)*(1-xi);
Dphi(2,1) = (1/4)*(1-eta);
Dphi(2,2) = (-1/4)*(1+xi);
Dphi(3,1) = (1/4)*(1+eta);
Dphi(3,2) = (1/4)*(1+xi);
Dphi(4,1) = (-1/4)*(1+eta);
Dphi(4,2) = (1/4)*(1-xi);


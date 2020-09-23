%---------------------------------------------------------------------------%                          
%                           Definition of the shape function derivatives    %
%                           in the master element.                          %
%---------------------------------------------------------------------------%

%Dphi_j is the derivative of the shape function  w.r.t. 
%the isoparametric coordinate.
%xi is the isoparametric coordinate
%nn is the number of nodes per element
function Dphi = shapefunctionFirstDerivative(xi,nn);

if( nn == 2)
    Dphi(1)= -1.0/2.0;  
    Dphi(2) = 1.0/2.0;
elseif( nn == 3)
    Dphi(1) = xi -0.5;
    Dphi(2) = -2*xi;
    Dphi(3) = xi +0.5;
elseif(nn == 4)
    Dphi(1) = (-9/16)*(3*xi^2-2*xi-(1/9));
    Dphi(2) = (27/16)*(3*xi^2-(2/3)*xi-1);
    Dphi(3) = (-27/16)*(3*xi^2+(2/3)*xi-1);
    Dphi(4) = (9/16)*(3*xi^2+(2*xi)-(1/9));
end


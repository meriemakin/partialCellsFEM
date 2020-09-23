%---------------------------------------------------------------------------%                          
%                           Definition of the shape functions               %
%                           in the master element.                          %
%---------------------------------------------------------------------------%

%phi are the shape functions 
%xi is the isoparametric coordinate
%nn is the number of nodes per element
function phi = shapefunction(xi,nn);

if(nn == 2)
    phi(1) = (1-xi)/2.0;
    phi(2) = (1+xi)/2.0;
elseif(nn == 3)
    phi(1) = 0.5*xi*(xi-1);
    phi(2) = -1.0*(xi+1)*(xi-1);
    phi(3) = 0.5*xi*(xi+1);
elseif(nn == 4)
    phi(1) = (-9/16)*(xi+(1/3))*(xi-(1/3))*(xi-1);
    phi(2) = (27/16)*(xi+1)*(xi-(1/3))*(xi-1);
    phi(3) = (-27/16)*(xi+1)*(xi+(1/3))*(xi-1);
    phi(4) = (9/16)*(xi+1)*(xi+(1/3))*(xi-(1/3));
end
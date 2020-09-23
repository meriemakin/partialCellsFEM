%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      %
%---------------------------------------------------------------------------%

%n is the total number of gauss points used for integration
%j is the jacobian
%A2 is a parameter of the DE.
%z is the material case
%nn is the number of nodes per element
function s = elementStiffnessCutCell(elem,A2,n,z,nn);

s = zeros(nn,nn);


%   USE OF ONE INTEGRATION POINT
    [x_g, w] = gaussianquadrature(n,1);
    jacob = jacobian2(elem,x_g);
    phi = shapefunction(x_g,nn);  
    Dphi = shapefunctionFirstDerivative(x_g,nn);
    [A1, dA1] = material(x_g, z);
    for l=1:nn
        for m=1:nn
            s(l,m) = s(l,m) + (dA1*Dphi(l)*phi(m)-A1*Dphi(l)*Dphi(m)*jacob+A2*phi(l)*phi(m)*(1/jacob))*w;
        end
    end


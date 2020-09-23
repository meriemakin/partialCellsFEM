%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      % 
%---------------------------------------------------------------------------%

%n is the total number of gauss points used for integration
%j is the jacobian
%A2 is a parameter of the DE.
%z is the material case
%nn is the number of nodes per element
function s = elementStiffness(elem,A1,n,z,nn);

s = zeros(nn,nn);
for i=1:n
    [x_g, w] = gaussianquadrature(n,i);
    jacob = jacobian2(elem,x_g);
    phi = shapefunction(x_g,nn);  
    Dphi = shapefunctionFirstDerivative(x_g,nn);
    for l=1:nn
        for m=1:nn
            s(l,m) = s(l,m) + (A1*Dphi(l)*Dphi(m)*jacob)*w;
        end
    end
end


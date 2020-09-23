%---------------------------------------------------------------------------%                          
%                           Definition of the element force vector          %
%---------------------------------------------------------------------------%

%n is the total number of gauss points used for integration
%j is the jacobian .
%K and L are factors in the DE.
%nn is the number of nodes per element
%elem is the considered element

function f = elementForceVector(n,K,L,nn,elem);

f = zeros(nn,1);
for i=1:n
    [x_g, w] = gaussianquadrature(n,i);
    phi = shapefunction(x_g,nn);
    jacob = jacobian2(elem,x_g);
    for l=1:nn
        f(l) = f(l) + phi(l)* K*w*(1/jacob);
    end
end
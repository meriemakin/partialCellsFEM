%---------------------------------------------------------------------------%                          
%                           Calculation of the energy norm                  % 
%---------------------------------------------------------------------------%

%n is the number of gauss points
%elem is the element for which the norm is calculated
%K  the factor in the DE
%L the length of the system
%A1 factor in the DE
%A2 factor in the DE
%bc1 left BC
%bc2 right BC
%u_n are the node displacements for the element elem
function norm = energyNormDifference(n,elem, K,L,A1,A2, bc1, bc2,u_n);

norm = 0.0;
for i=1:n
    [x_g, w] = gaussianquadrature(n,i);
    phi = shapefunction(x_g,length(elem));
    Dphi = shapefunctionFirstDerivative(x_g,length(elem));
    x = mapping(x_g,phi,elem); 
    jacob = jacobian2(elem,x_g);
    derivativeExact = derivativeExactSolution(x, K,L,A1,A2, bc1, bc2);
    derivativeFem = derivativeFEMSolution(Dphi,u_n);
    norm = norm + derivativeFem * A1 * derivativeFem * jacob *w - ...
    2* derivativeFem * A1 * derivativeExact * w + ...
    derivativeExact * A1* derivativeExact * (1/jacob) * w;
end






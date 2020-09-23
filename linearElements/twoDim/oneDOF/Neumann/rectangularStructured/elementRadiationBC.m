%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%n is the total number of gauss points used for integration
%nn is the number of nodes per element
function s = elementRadiationBC(elem,n,nn,A,h);

s = zeros(nn,nn);

for j=1:n
    [xi_g,eta_g,w1,w2] = gaussianquadrature(n,1,j);
    phi = shapefunction(1.0, eta_g);  
    Dphi = shapefunctionFirstDerivative(1.0,eta_g);
    F = deformationGradient(elem,Dphi);
    jacob = jacobian(F);
    s = s + A*h*phi'*phi*jacob*w2;        
end



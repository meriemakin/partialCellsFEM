%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%n is the total number of gauss points used for integration
%nn is the number of nodes per element
function s = elementStiffness(elem,n,nn);

s = zeros(nn,nn);

for i=1:n
  for j=1:n
    [xi_g,eta_g,w1,w2] = gaussianquadrature(n,i,j);
    phi = shapefunction(xi_g, eta_g);  
    Dphi = shapefunctionFirstDerivative(xi_g,eta_g);
    F = deformationGradient(elem,Dphi);
    jacob = jacobian(F);
    s = s + Dphi*((F'*F)^-1)*Dphi'*jacob*w1*w2;        
  end
end



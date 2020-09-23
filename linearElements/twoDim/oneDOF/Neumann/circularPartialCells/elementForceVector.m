%---------------------------------------------------------------------------%                          
%                           Definition of the element force vector          %
%---------------------------------------------------------------------------%

%n is the total number of gauss points used for integration
%nn is the number of nodes per element
%elem is the considered element
%flags are the bc flags
%pen is the penalty parameter
%T0 temperature boundary condition
function f = elementForceVector(n,nn,elem,A,rhs);

f = zeros(nn,1);

for i=1:n
  for j=1:n
    [xi_g,eta_g,w1,w2] = gaussianquadrature(n,i,j);
    phi = shapefunction(xi_g, eta_g); 
    Dphi = shapefunctionFirstDerivative(xi_g,eta_g);
    F = deformationGradient(elem,Dphi);
    jacob = jacobian(F);
    f = f + (rhs/A) * phi'*jacob*w1*w2;
  end
end


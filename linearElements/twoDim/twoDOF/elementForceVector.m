%---------------------------------------------------------------------------%                          
%                           Definition of the element force vector          %
%---------------------------------------------------------------------------%

%n is the total number of gauss points used for integration
%nn is the number of nodes per element
%elem is the considered element
function f = elementForceVector(n,nn,elem,rhs,flag);

f = zeros(nn*2,1);

if(flag == 1 || flag == 2)
    for i=1:n
       for j=1:n
       [xi_g,eta_g,w1,w2] = gaussianquadrature(n,i,j);
        phi = shapefunction(xi_g, eta_g); 
        N = shapefunctionMatrix(phi);
        Dphi = shapefunctionFirstDerivative(xi_g,eta_g);
        F = deformationGradient(elem,Dphi);
        jacob = det(F);
        f = f + N'*rhs*jacob*w1*w2;
        end
    end
end

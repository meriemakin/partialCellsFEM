%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%n is the total number of gauss points used for integration
%nn is the number of nodes per element
%A is coefficient in differential equation to be solved
function s = elementStiffness(elem,n,nn,A,flag);

s = zeros(nn*2,nn*2);

if(flag == 1 || flag == 2)
	for i=1:n
  		for j=1:n
    			[xi_g,eta_g,w1,w2] = gaussianquadrature(n,i,j);
    			phi = shapefunction(xi_g, eta_g);  
    			Dphi = shapefunctionFirstDerivative(xi_g,eta_g);
    			F = deformationGradient(elem,Dphi);
    			DphiR = shapefunctionFirstDerivativeReal(Dphi, F);
    			bm =  BMatrix(DphiR);
    			jacob = det(F);
    			s = s + A*bm'*bm*jacob*w1*w2;        
  		end
	end
end


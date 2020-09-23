%---------------------------------------------------------------------------%                          
%                           Derivative of the FEM solution                  %
%---------------------------------------------------------------------------%

%n is the number of nodes per element
%u_n 
function du_n = derivativeFEMSolution(n,u_n, elementL);

nn = length(u_n);

du_n = zeros(nn,1);

%this is done only for element with two nodes

jacobian = 2.0/ elementL;

for j = 1:(nn - 1)
	for i=1:n
             Dphi = shapefunctionFirstDerivative(-1.0,n);
    	     du_n(j,1) = du_n(j,1) + u_n(j+(i-1),1)*Dphi(i)*jacobian;
	end
end

du_n(nn,1) = du_n(nn-1,1);


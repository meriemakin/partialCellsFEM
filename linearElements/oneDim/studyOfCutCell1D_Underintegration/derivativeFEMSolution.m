%---------------------------------------------------------------------------%                          
%                           Derivative of the FEM solution                  %
%---------------------------------------------------------------------------%

%Dphi are the shape functions derivative
%u_n 
function du = derivativeFEMSolution(Dphi,u_n);

du = 0.0;

for i=1:length(Dphi)
    du = du + u_n(i)*Dphi(i);
end


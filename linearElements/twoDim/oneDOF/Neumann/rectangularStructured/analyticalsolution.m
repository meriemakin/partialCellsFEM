%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%analytical solution for laplace(u) = rhs calculated using a green's function with homogeneoous zero boundary conditions.

function [anasol] = analyticalsolution(y_n,bc);


b = max(y_n);

anasol = zeros(length(y_n),1);

for i=1:length(y_n)
      y = y_n(i);
      anasol(i,1) =  ((bc*y)/b);
end


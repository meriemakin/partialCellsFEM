%------------------------------------------------------------------------
%                   Calculation of the error made by the algorithm 
%------------------------------------------------------------------------

%u_n is the numerical solution
%u_nE is the analytical solution
function [nin,err] = errorNorm(u_n, u_nE);

err = 0.0;

nin = length(u_n) - 2;

for i=2:(nin+1)

     err = err + (abs(u_n(i)-u_nE(i)))^2;

end

%err = sqrt(err);
err = err/nin;
err = sqrt(err);

%---------------------------------------------------------------------------%                          
%                         Analytical Solution to A Circle With Dirichlet BC %
%---------------------------------------------------------------------------%

%analytical solution for laplace(u) = rhs calculated using a green's function with homogeneoous zero boundary conditions.

function [anasol] = analyticalsolution(x_n,y_n,f,A,n);


a = max(x_n);
b = max(y_n);

anasol = zeros(length(x_n),1);

for i=1:length(x_n)
      x = x_n(i);
      y = y_n(i);
          
     for m=1:n
         for k=1:n
         anasol(i,1) = anasol(i,1) +  ((1-cos(k*pi))*(1-cos(m*pi))*(sin((k*pi*x)/a))*(sin((m*pi*y)/b)))/(m*k*((m^2/b^2)+(k^2/a^2)));
         end
     end
end

anasol = ((4*f)/(A*pi^4)) * anasol;

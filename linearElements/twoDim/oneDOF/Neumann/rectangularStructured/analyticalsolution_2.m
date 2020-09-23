%---------------------------------------------------------------------------%                         
%                         Analytical Solution to A Circle With Dirichlet BC % 
%---------------------------------------------------------------------------%

%analytical solution for laplace(u) = rhs calculated using a green's function with homogeneoous zero boundary conditions.

function [anasol] = analyticalsolution_2(x_n,y_n,bc,h,n);


a = max(x_n);
b = max(y_n);

anasol = zeros(length(x_n),1);
%for h*Lx = 0.2
%alphas are solutions of alpha*tan(alpha) = h*Lx)
%solve alpha*tan(alpha) = h*Lx, numerically:
f1 = @(x) (x*tan(x)-(h*a));
alphas = [];
for i=0.1:n
	pt = fzero(@(x) f1(x),i*pi);  % The x-coord
	alphas = [alphas;pt];
end
alphas'
alphas = alphas/a;

for i=1:length(x_n)
      x = x_n(i);
      y = y_n(i);
      for m=1:n
	 al = alphas(m);
	 fraction = cosh(al*(b-y))/cosh(al*b);
         anasol(i,1) = anasol(i,1) +  fraction*(cos(al*x))/(((al^2+h^2)*a+h)*cos(al*a));
     end
end

anasol = 2*h*bc * anasol;

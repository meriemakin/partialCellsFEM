%---------------------------------------------------------------------------%                          
%                           Definition of the element force vector          %
%---------------------------------------------------------------------------%

%n is the total number of gauss points used for integration
%j is the jacobian .
%K is the factor in the DE.
%L is the real end of the domain
%nn is the number of nodes per element
%elem is the considered element
%bc2 is the prescribed boundary condition at the right end of the domain

function f = elementForceVectorCutCell(n,L,bc2,nn,elem,penalty);

f = zeros(nn,1);
X1 = elem(1);
X2 = elem(length(elem));
A = [X1, 1; X2, 1];
R = [-1 1];
sol = A/R;
%a = -1/((3*0.5*X1)+(0.5*X2))
%b = a*0.5*(X1 + X2)
x_g = (sol(1) * L) + sol(2);
phi = shapefunction(x_g,nn);
jacob = jacobian(elem);
for l=1:nn
    %f(l) = f(l) + 1000*phi(l)*(1/jacob)*bc2;
    %f(l) = f(l) + phi(l)*(1/jacob)*bc2;
     

     %f(l) = f(l) +  penalty*phi(l)*bc2;
     %In order to have positive entries on the diagonal
     %implemented is - ... instead of + , but it leads
     %to the same results. Adding a zero or substracting a 
     %zero is the same at the end.
     f(l) = f(l) -  penalty*phi(l)*bc2;
end




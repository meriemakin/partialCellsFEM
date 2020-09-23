%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      %
%---------------------------------------------------------------------------%

%n is the total number of gauss points used for integration
%j is the jacobian
%A2 is a parameter of the DE.
%z is the material case
%nn is the number of nodes per element
%L is the real end of the domain
function s = elementStiffnessCutCell(elem,n,nn,L,penalty);

s = zeros(nn,nn);

jacob = jacobian(elem);     
%evaluate the shape functions at the real end of the domain
X1 = elem(1);
X2 = elem(length(elem));
A = [X1, 1; X2, 1];
R = [-1 1];
sol = A/R;
%a = -1/((3*0.5*X1)+(0.5*X2))
%b = a*0.5*(X1 + X2)
x_g = (sol(1) * L) + sol(2);
phi = shapefunction(x_g,nn)
for l=1:nn
    for m=1:nn
        %s(l,m) = s(l,m) + 1000*phi(l)*phi(m)*(1/jacob);
        %s(l,m) = s(l,m) + phi(l)*phi(m)*(1/jacob);
        

        %s(l,m) = s(l,m) + penalty*phi(l)*phi(m);
        %In order to have positive entries on the diagonal
        %implemented is - ... instead of + , but it leads
        %to the same results. Adding a zero or substracting a 
        %zero is the same at the end.
        s(l,m) = s(l,m) - penalty*phi(l)*phi(m);
    end
end


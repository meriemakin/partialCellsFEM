%---------------------------------------------------------------------------%                          
%                           Comparing of error to tolerance                 %
%---------------------------------------------------------------------------%

%x_n are the nodes
%nn is the number of nodes per element
%tol is the error tolerance
%n is the number of gauss points
%K  the factor in the DE
%L the length of the system
%A1 factor in the DE
%A2 factor in the DE
%bc1 left BC
%bc2 right BC
%u_n are the fem solution on nodes
function errors = errorCalc(x_n,nn,tol,n,K,L,A1,A2,bc1,bc2, u_n);

elem = elements(x_n,nn);
[nelements, nnodes] = size(elem);
errors = zeros(nelements,1);
error1 = 0.0;
error2 = 0.0;

for i= 1: nelements
    element = elem(i,:);
    error1 = error1 + energyNormExactFunction(n,element, K,L,A1,A2, bc1, bc2)
end

for i=1:nelements
    element = elem(i,:);
    u_ne = zeros(nn,1);
    for j = 1:nn
        u_ne(j,1) = u_n(j+(i-1)*(nn-1));
    end
    error2 = error2 + energyNormDifference(n,element, K,L,A1,A2, bc1, bc2,u_ne);
end

errors = sqrt(error2/error1);

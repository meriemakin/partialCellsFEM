%---------------------------------------------------------------------------%                          
%                           generating nodes                                %
%---------------------------------------------------------------------------%

%L is the length of the system
%equidistant is true if we want equal sized elements, and false if not.
%n is the total number of elements
%nn is the number of nodes per element
function x_n = vertices(L,equidistant,n,nn);

%total number of nodes
tn = n*(nn-1)+1;
x_n = zeros(tn,1);
for i=1:length(x_n)
    x_n(i) = (L/(tn-1))*(i-1);
end



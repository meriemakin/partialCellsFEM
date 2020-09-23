%---------------------------------------------------------------------------%                          
%                           generating nodes                                %
%---------------------------------------------------------------------------%

%xs is the starting coordinate of the domain
%L is the length of the system
%n is the total number of elements
%nn is the number of nodes per element
function x_n = nodes_vers2(xs,xe,n,nn);

%total number of nodes
tn = n*(nn-1)+1;
x_n = zeros(tn,1);
L = xe - xs;
for i=1:length(x_n)
    x_n(i) = xs + (L/(tn-1))*(i-1);
end



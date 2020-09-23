%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%x_n is the nodes coordinates
%n number of nodes per element
function elem = elements(x_n, n);

l = (length(x_n)-1)/(n-1);
elem = zeros(l,n);

for i=1:l
    for j=1:n
        elem(i,j) = x_n((n-1)*(i-1)+j);
    end
end

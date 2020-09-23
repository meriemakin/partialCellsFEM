%---------------------------------------------------------------------------%                          
%                           generating centrois                             %
%---------------------------------------------------------------------------%
%x_n are the vertices coordinates
function c_n = centroids(x_n);

n = length(x_n);
m = n -1;

c_n = zeros(m,1);
for i=2:n
    c_n(i-1,1)=(x_n(i-1,1)+x_n(i,1))*0.5;
end

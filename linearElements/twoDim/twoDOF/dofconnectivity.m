%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%Nr is the total number of elements in r-direction
%Ntheta is the total number of element in theta-direction
%nxis the number of nodes per elements in x-direction
%ny is the number of nodes per elements in y-direction
function dofconn = dofconnectivity(conn);

[r,c] = size(conn);
dofconn = zeros(r , c * 2);

for i=1:r
    for j=1:c
        dofconn(i,(j-1)*2 + 1) = (conn(i,j) -1 )*2 + 1; 
        dofconn(i,(j-1)*2 + 2) = (conn(i,j) -1 )*2 + 2;
    end
end

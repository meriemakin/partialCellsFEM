%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%Nr is the total number of elements in r-direction
%Ntheta is the total number of element in theta-direction
%nxis the number of nodes per elements in x-direction
%ny is the number of nodes per elements in y-direction
function conn = connectivity(Nr, Ntheta,ny,nx);

conn = zeros(Nr*Ntheta, ny*nx);

for k=1:Ntheta
     for l=1:Nr
         conn(Nr*(k-1)+l,1:4) = [(Nr+1)*(k-1)+l,(Nr+1)*(k)+l,(Nr+1)*(k)+(l+1), (Nr+1)*(k-1)+(l+1)];
         %conn(Nr*(k-1)+l,1:4) = [(Nr+1)*(k-1)+l,(Nr+1)*(k-1)+(l+1), (Nr+1)*(k)+(l+1),(Nr+1)*(k)+(l)];
     end
end

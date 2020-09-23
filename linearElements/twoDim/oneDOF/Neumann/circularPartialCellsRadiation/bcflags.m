%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%Nr is the total number of elements in r-direction
%Ntheta is the total number of element in theta-direction
%nxis the number of nodes per elements in x-direction
%ny is the number of nodes per elements in y-direction
function flags = bcflags(elems);

[n,r,c] = size(elems);
flags = zeros(c,r);
%0 no BCs
%1 temperature BCs
%2 Flux BCs
for i=1:c
    element = elems(:,:,i);
    for j=1:r
        x1 = element(1,j);
        y1 = element(2,j);
        if(j < r)
            x2 = element(1,j+1);
            y2 = element(2,j+1);
        elseif ( j == r)
            x2 = element(1,1);
            y2 = element(2,1);
        end
        if(y1 == 0 && y2 == 0 && x1>0 && x2>0)
            flags(i,j) = 2;
        elseif (((y1 - 0) < 10^(-10)) && ((y2 - 0) < 10^(-10)) && x1<0 && x2<0)
            flags(i,j) = 1;
        end
    end
end

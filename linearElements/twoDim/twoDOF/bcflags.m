%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%elems are the elements
%elemflags: 0 outside of domain, 1 totally inside of domain, 2 partially inside of domain
%xC and yC are the coordinates of the center of the circular domoain
%radius is the radius of the circle
function flags = bcflags(elems,elemFlags,xC,yC, radius);

[n,r,c] = size(elems);
flags = zeros(c,1);

%0 no BCs
%1 BCs

%Canonical problem
%The right upper quarter of the circular physical domain is subject to boundary conditions

for i=1:c
    element = elems(:,:,i);
    flagE = elemFlags(i,1);
    count = 0;
    if(flagE == 2)
       for j=1:r
           x = element(1,j);
           y = element(2,j);
           if(x > xC && y > yC )
              count = count + 1;
           end
        end
        if(count >=1 ) 
           flags(i,1) = 1;
        end
    end
end

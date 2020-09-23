%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%Nr is the total number of elements in r-direction
%Ntheta is the total number of element in theta-direction
%nxis the number of nodes per elements in x-direction
%ny is the number of nodes per elements in y-direction
function numelem = radiationflags(elems,conn,xmax);


numelem = [];
[n,r,c] = size(elems);

for i=1:c

    element  = elems(:,:,i);
    connElem = conn(i,:);

    for j=1:r
        x = element(1,j);
        y = element(2,j);
        if(x == xmax)
           numelem = [numelem;i];
	   break;
        end
    end
end

numelem = sort(unique(numelem));

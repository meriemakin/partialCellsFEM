%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%Nr is the total number of elements in r-direction
%Ntheta is the total number of element in theta-direction
%nxis the number of nodes per elements in x-direction
%ny is the number of nodes per elements in y-direction
function dofs = bcflags(elems,conn,xmin,xmax,ymin,ymax);


dofs = [];
[n,r,c] = size(elems);
for i=1:c

    element  = elems(:,:,i);
    connElem = conn(i,:);

    for j=1:r
        x = element(1,j);
        y = element(2,j);
        if(x == xmin || x == xmax || y == ymin || y == ymax)
           dofs = [dofs;connElem(j)];
        end
    end
end

dofs = sort(unique(dofs));

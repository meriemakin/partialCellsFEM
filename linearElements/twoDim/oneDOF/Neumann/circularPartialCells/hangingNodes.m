%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%Nr is the total number of elements in r-direction
%Ntheta is the total number of element in theta-direction
%nxis the number of nodes per elements in x-direction
%ny is the number of nodes per elements in y-direction
function hN = hangingNodes(elem,xC,yC,radius,connElem);


hN = [];
[r,c] = size(elem);


for j=1:c
    x = elem(1,j);
    y = elem(2,j);
    if(sqrt((xC- x)^2.0 + (yC - y)^2.0) > radius)
       hN=[hN;connElem(j),x, y];
    elseif(sqrt((xC- x)^2.0 + (yC - y)^2.0) == radius)
       hN=[hN;connElem(j),x, y];
    end
end



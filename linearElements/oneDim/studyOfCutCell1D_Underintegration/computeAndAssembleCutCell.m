%---------------------------------------------------------------------------%                          
%                           Assembling of the element stiffness matrix      %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

%x_n the nodes positions
%A2 is a coefficient in the DE
%z is the type of the material behavior (constant or nonconstant)
%n is the total number of gauss points per element
%K is the factor in the DE
%L is the system length
%nn is the number of nodes per element
function [s,f] = computeAndAssembleCutCell(x_n,A2,z, n,K,L, nn);

elem = elements(x_n,nn)
s = zeros(length(x_n), length(x_n));
f = zeros(length(x_n),1);
[nelements, numberofnodes] = size(elem);
for i =1:nelements
    if(i ~= nelements)
      sElem = elementStiffness(elem(i,:),A2,n,z,nn);
    else
      sElem = elementStiffnessCutCell(elem(i,:),A2,n,z,nn);  
    end
    
    for j = 1:nn
        for k = 1:nn
            s((nn-1)*(i-1)+j,(nn-1)*(i-1)+k) = s((nn-1)*(i-1)+j,(nn-1)*(i-1)+k) + sElem(j,k);
        end
    end
    if(i ~= nelements)
      fElem = elementForceVector(n,K,L,nn,elem(i,:));
    else
      fElem = elementForceVectorCutCell(n,K,L,nn,elem(i,:));  
    end    
    for j = 1:nn
        f((nn-1)*(i-1)+j,1) = f((nn-1)*(i-1)+j,1) + fElem(j,1);
    end
end

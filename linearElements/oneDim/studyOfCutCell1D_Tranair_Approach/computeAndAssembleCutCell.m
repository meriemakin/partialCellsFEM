%---------------------------------------------------------------------------%                          
%                           Assembling of the element stiffness matrix      %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

%x_n the nodes positions
%A2 is a coefficient in the DE
%z is the type of the material behavior (constant or nonconstant)
%n is the total number of gauss points per element
%K is the factor in the DE
%L is real end of the system
%nn is the number of nodes per element
function [s,f] = computeAndAssembleCutCell(x_n,A2,z, n,K,L, nn,bc2,p,kappa);

elem = elements(x_n,nn);
s = zeros(length(x_n), length(x_n));
f = zeros(length(x_n),1);
[nelements, numberofnodes] = size(elem);
for i =1:nelements    
    sElem = elementStiffness(elem(i,:),A2,n,z,nn);
    %This has been added to appropriately consider the contribution of the partial differetntial equation to the real domain. Since the finite element domain has been extended to match the finite volume unknowns positions, only one half of the first finite element corresponds to the real physical domain. 
    if(i == 1)
      sElem = 0.5*sElem;
    end
 
    if(i == nelements)
      temp = elementStiffnessCutCell(elem(i,:),n,nn,L,p);
      sElem = kappa*sElem;
      sElem = sElem + temp;
    end
    
    for j = 1:nn
        for k = 1:nn
            s((nn-1)*(i-1)+j,(nn-1)*(i-1)+k) = s((nn-1)*(i-1)+j,(nn-1)*(i-1)+k) + sElem(j,k);
        end
    end
    
    fElem = elementForceVector(n,K,nn,elem(i,:));    
    
    if(i == 1)
      fElem = 0.5*fElem;
    end

    if(i == nelements)
      temp = elementForceVectorCutCell(n,L,bc2,nn,elem(i,:),p);
      fElem = kappa*fElem;
      fElem = fElem+ temp; 
    end    
    
    for j = 1:nn
        f((nn-1)*(i-1)+j,1) = f((nn-1)*(i-1)+j,1) + fElem(j,1);
    end
end

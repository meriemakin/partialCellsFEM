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
%b is the element that one wants to begin with solving the PDE
%it is 1 if we would like to proceed with the first element and 2 if we want to
%not use the first element but impose a constraint instead

function [s,f] = computeAndAssembleCutCell(x_n,A1,z, n,K,L, nn,b,alpha,flag);

elem = elements(x_n,nn);
s = zeros(length(x_n), length(x_n));
f = zeros(length(x_n),1);
[nelements, numberofnodes] = size(elem);
for i = b:(nelements)
    sElem = elementStiffness(elem(i,:),A1,n,z,nn);
    if (flag == 1)
    	if (i == b)
      	 	sElem = 0.5* sElem; 
    	end
    	if (i == nelements)
        	sElem = alpha *  sElem;
    	end 
    end  
    for j = 1:nn
        for k = 1:nn
            s((nn-1)*(i-1)+j,(nn-1)*(i-1)+k) = s((nn-1)*(i-1)+j,(nn-1)*(i-1)+k) + sElem(j,k);
        end
    end

    fElem = elementForceVector(n,K,L,nn,elem(i,:));
    if (flag == 1)
    	if (i == b)
       		fElem = 0.5 * fElem;
    	end
    	if (i == nelements)
       		fElem = alpha * fElem;
    	end
    end
    for j = 1:nn
        f((nn-1)*(i-1)+j,1) = f((nn-1)*(i-1)+j,1) + fElem(j,1);
    end
end




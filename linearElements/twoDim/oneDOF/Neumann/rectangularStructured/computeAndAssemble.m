%---------------------------------------------------------------------------%                          
%                           Assembling of the element stiffness matrix      %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

%elem are the elements
%conn is the connectivity
%ng is the number of gauss points in one direction
%rhs the value of the function f in the PDE
%tnn the total number of nodes in the domain
function [solution] = computeAndAssemble(elem,conn,ng,A,rhs,tnn,bc,x_n,y_n,h);

[nelements, numberofnodes] = size(conn);
%s = zeros(tnn,tnn);
s = sparse(tnn, tnn);
f = zeros(tnn,1);

%numelem = radiationflags(elem,conn,max(x_n));

for i = 1:nelements 

    	sElem = elementStiffness(elem(:,:,i),ng,numberofnodes,A);
   
    	connElem = conn(i,:);
    
    	for j = 1:numberofnodes
        	for k = 1:numberofnodes
            		s(connElem(j),connElem(k)) = s(connElem(j),connElem(k)) + sElem(j,k);
        	end
    	end
   	
	%found = searchindex(numelem,i);
	%if(found == 1)
	%	sradiation = elementRadiationBC(elem(:,:,i),ng,numberofnodes,A,h);
	%	
	%	for j = 1:numberofnodes
	%                for k = 1:numberofnodes
        
        %	                s(connElem(j),connElem(k)) = s(connElem(j),connElem(k)) + sradiation(j,k);
        %       	end
        %	end

	%end 
    	
	fElem = elementForceVector(ng,numberofnodes,elem(:,:,i),rhs);

    	for j = 1:numberofnodes
        	f(connElem(j),1) = f(connElem(j),1) + fElem(j,1);
    	end
end

%apply boundary conditions to the nodes that lie on the boundary


%apply dirichlet boundary conditions
dofs_1 = bcflags(elem,conn,min(y_n));
dofs_2 = bcflags(elem,conn,max(y_n));
[s,f] = applyBC(s,f,bc,dofs_1,dofs_2);
dofs = [dofs_1;dofs_2];
r1 = length(dofs_1);
solution = s\f;

for i=1:length(dofs)
    if(i>r1)	
    	solution = [solution(1:dofs(i)-1);bc;solution(dofs(i):end)];
    else
	solution = [solution(1:dofs(i)-1);0;solution(dofs(i):end)];	
    end
end


dofs_2 = bcflags(elem,conn,max(y_n));


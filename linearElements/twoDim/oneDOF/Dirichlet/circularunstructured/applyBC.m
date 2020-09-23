%---------------------------------------------------------------------------%                          
%                           Apply the bc on the assembled stiffness matrix  %
%                           and force vector                                %
%                           and the element force vector                    % 
%---------------------------------------------------------------------------%

%s1 is the assembled stiffness matrix
%f1 is the assembled force vector
%bc2 is the right bc
%bc1 is the left bc
function [s,f] = applyBC(s,f,bc,dofs);

dofs_copy = dofs;
r = length(dofs_copy);

for i =1:r
        s(dofs_copy(i),:) = [];
	f(dofs_copy(i)) = [];
	for j=1:length(f)
    		f(j,1)= f(j,1) - s(j,dofs_copy(i))*bc;
	end
	s(:,dofs_copy(i)) =[];
        dofs_copy = dofs_copy - ones(1,length(dofs_copy));
end

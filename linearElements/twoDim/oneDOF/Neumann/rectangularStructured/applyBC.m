%---------------------------------------------------------------------------%                          
%                           Apply the bc on the assembled stiffness matrix  %
%                           and force vector                                %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

%s1 is the assembled stiffness matrix
%f1 is the assembled force vector
%bc2 is the right bc
%bc1 is the left bc
function [s,f] = applyBC(s,f,bc,dofs_1,dofs_2);

r1 = length(dofs_1);
dofs_copy = [dofs_1;dofs_2];
r = length(dofs_copy);

for i =1:r
        s(dofs_copy(i),:) = [];
	f(dofs_copy(i)) = [];
	for j=1:length(f)
		if(i >r1)
			f(j,1)= f(j,1) - s(j,dofs_copy(i))*bc;
		end	
	end
	s(:,dofs_copy(i)) =[];
        dofs_copy = dofs_copy - ones(length(dofs_copy),1);
end

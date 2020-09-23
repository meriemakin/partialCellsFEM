%s1 is the assembled stiffness matrix
%f1 is the assembled force vector
%bc2 is the right bc
%bc1 is the left bc
function [s,f] = applyBCDirichlet(s,f,bc,dof);

s(dof,:) = [];
f(dof) =[];

for i=1:length(f)
    f(i,1)= f(i,1) - s(i,dof)*bc;
end

s(:,dof) =[];


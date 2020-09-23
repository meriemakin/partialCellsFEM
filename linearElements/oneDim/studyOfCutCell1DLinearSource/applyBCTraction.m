%---------------------------------------------------------------------------%                          
%                           Apply the bc on the assembled stiffness matrix  %
%                           and force vector                                %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

%s1 is the assembled stiffness matrix
%f1 is the assembled force vector
%bc2 is the right bc
%bc1 is the left bc
function [s,f] = applyBCTraction(s,f,bc1, bc2);

%delete the first and last row
s(1,:) = [];
s(length(s)-1,:) = [];
f(1) =[];
f(length(f)) = [];

for i=1:length(f)
    f(i,1)= f(i,1) - s(i,1)*bc1 - s(i,length(s))*bc2;
end

%delete the first and last column
s(:,1) =[];
s(:,length(s)) = [];


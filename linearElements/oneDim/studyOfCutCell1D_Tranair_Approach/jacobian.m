%---------------------------------------------------------------------------%                          
%                           Definition of the jacobian                      %
%---------------------------------------------------------------------------%

%elem is the element for which the jacobian is computed
function J = jacobian(elem);

J = 2/(elem(length(elem))-elem(1));

if(J < 0)
    disp('WARNING!!!!!!!Jacobian is negative!!!!!!');
end
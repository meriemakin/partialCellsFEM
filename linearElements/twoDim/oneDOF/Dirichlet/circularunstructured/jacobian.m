%---------------------------------------------------------------------------%                          
%                           Definition of the jacobian                      %
%---------------------------------------------------------------------------%

%xi and eta are the coordinates of the gauss point
%elem is the element for which the jacobian is computed
function J = jacobian(F);

J = det(F);
if(J <= 0)
disp('Jacobian is negative!');
end

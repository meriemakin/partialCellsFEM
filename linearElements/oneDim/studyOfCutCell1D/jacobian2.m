%---------------------------------------------------------------------------%                          
%                           Definition of the jacobian                      %
%---------------------------------------------------------------------------%

%elem is the element for which the jacobian is computed
%x_g is the coordinate of the gauss point
function J = jacobian2(elem,x_g);

%nn = length(elem);
%J_old = 2/(elem(length(elem))-elem(1))
%Dphi = shapefunctionFirstDerivative(x_g,nn);
%J = 1/(elem*Dphi')

J = 2/(elem(length(elem))-elem(1));
if(J < 0)
    disp('WARNING!!!!!!!Jacobian is negative!!!!!!');
end

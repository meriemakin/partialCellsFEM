%---------------------------------------------------------------------------%                          
%                           Definition of the jacobian                      %
%---------------------------------------------------------------------------%

%xi and eta are the coordinates of the gauss point
%elem is the element for which the jacobian is computed
% in elem are the nodes coordinates stored
%Dphi are the shape functions derivatives
function F = deformationGradient(elem,Dphi);

F = elem*Dphi;

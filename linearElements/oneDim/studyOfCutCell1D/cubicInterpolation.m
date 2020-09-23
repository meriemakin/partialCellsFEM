%---------------------------------------------------------------------------%                          
%                           Apply the bc on the assembled stiffness matrix  %
%                           and force vector                                %
%                           and the element force vector                    % 
%---------------------------------------------------------------------------%

function [weight1, weight2, weight3, weight4,rhs] = cubicInterpolation(x1,x2,x3,x4,u4,x5);

Deltax = x5 - x3;
rhs = u4;
weight1 = - ((x4-x2)*(x4-x3)*(x4-x5))/(6*Deltax^3);
weight2 = (x4*(x4-x3)*(x4-x5))/(2*Deltax^3);
weight3 = -(x4*(x4-x2)*(x4-x5))/(2*Deltax^3);
weight4 = (x4*(x4-x2)*(x4-x3))/(6*Deltax^3);

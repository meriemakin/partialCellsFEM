%---------------------------------------------------------------------------%                          
%                           Apply the bc on the assembled stiffness matrix  %
%                           and force vector                                %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

function [weight1, weight2, weight4,rhs] = quadInterpolationDerivative(x1,x2,x3,du3);

deltaX = x2 - x1;

weight1 = (2 * (x3 - x1) - 3 * deltaX) / (2* deltaX^2.0);
weight2 = -1.0 * (2 * ( x3 - x1 ) - 2 * deltaX) / (deltaX^2.0);
weight4 =  (2 * (x3 - x1) - deltaX) / (2 * deltaX^2.0);

rhs = du3;


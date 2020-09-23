%---------------------------------------------------------------------------%                          
%                           Apply the bc on the assembled stiffness matrix  %
%                           and force vector                                %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

function [weight1, weight2, rhs] = linearInterpolationDerivative(x1,x2,du3);

deltaX = x2 - x1;

weight1 = -1/deltaX;
weight2 =  1/deltaX;

rhs = du3;


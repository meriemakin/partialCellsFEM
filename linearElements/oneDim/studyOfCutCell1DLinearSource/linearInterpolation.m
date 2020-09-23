%---------------------------------------------------------------------------%                          
%                           Linear Interpolation to impose the BC           %
%---------------------------------------------------------------------------%

function [weight1, weight2, rhs] = linearInterpolation(x1,x2,x3,u3);


Deltax = x2 -x1;
newx = x3 - x1;


weight1 = -1.0 * (newx - Deltax)/Deltax;
weight2 = newx / Deltax;

rhs = u3;


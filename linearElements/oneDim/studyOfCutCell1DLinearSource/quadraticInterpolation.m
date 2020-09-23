%---------------------------------------------------------------------------%                          
%                           Apply the bc on the assembled stiffness matrix  %
%                           and force vector                                %
%                           and the element force vector                    % 
%---------------------------------------------------------------------------%

function [weight1, weight2, weight4,rhs] = quadraticInterpolation(x1,x2,x3,u3,x4);

alpha_a1 = (1/(x3^2-x1*x3-x2*x3+x1*x2))*(((x3-x1)/(x2-x1))-1);
alpha_a2 = (-1/(x3^2-x1*x3-x2*x3+x1*x2))*((x3-x1)/(x2-x1));
c_a = u3*(1/(x3^2-x1*x3-x2*x3+x1*x2));

alpha_b1 = -((alpha_a1*(x1+x2)) + (1/(x2-x1)));
alpha_b2 = ( (1/(x2-x1)) - alpha_a2*(x1+x2));
c_b = -1*c_a*(x1+x2);

alpha_c1 = 1 - alpha_a1*x1^2 - alpha_b1*x1;
alpha_c2 = -1.0*alpha_a2*x1^2 - alpha_b2*x1;
c_c = -1.0*c_a*x1^2-c_b*x1;

weight1 = -1.0*(alpha_a1*x4^2 + alpha_b1 * x4 + alpha_c1);
weight2 = -1.0*(alpha_a2*x4^2 + alpha_b2 * x4 + alpha_c2);
weight4 = 1;

rhs = c_a*x4^2 + c_b*x4 + c_c;


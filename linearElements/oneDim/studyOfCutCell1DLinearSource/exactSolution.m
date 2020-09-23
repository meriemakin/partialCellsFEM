%---------------------------------------------------------------------------%                          
%                           Analytical solution                             % 
%---------------------------------------------------------------------------%

%x_n the nodes positions
%K  the factor in the DE
%L the length of the system
%A1 factor in the DE
%A2 factor in the DE
%bc1 left BC
%bc2 right BC
function u = exactSolution(x_n, k,L,A1,A2, bc1, bc2);

u = 0.0;


c2 = bc1*A1;
c1 = (bc2- ((k/(2*A1))*L^2) - (c2/A1))*(A1/L);
u = (k/(2*A1))*x_n.^2+(c1/A1)*x_n+(c2/A1);

%u = 50*x_n.^2 - 49*x_n;
%---------------------------------------------------------------------------%                          
%                           Derivative of analytical solution               %
%---------------------------------------------------------------------------%

%x_n the nodes positions
%K  the factor in the DE
%L the length of the system
%A1 factor in the DE
%A2 factor in the DE
%bc1 left BC
%bc2 right BC
function du = derivativeExactSolution(x_n, k,L,A1,A2, bc1, bc2);

du = 0.0;

%homogeneous part of the solution
if(A2 ~= 0 && k~=0)
    r1 = sqrt(abs(A2/A1));
    r2 = -sqrt(abs(A2/A1));
    if( A2 * A1 < 0)
        c2 = (bc2- bc1*exp(r1*L))/(exp(r2*L)-exp(r1*L));
        c1 = bc1 - c2; 
        du = du+ c1 * r1* exp(r1*x_n) + c2 * r2 * exp(r2*x_n);
    elseif( A2 * A1 > 0)
        c2 = (bc2- bc1*cos(r1*L))/sin(r1*L);
        c1 = bc1; 
        du = du - c1 * r1 * sin(r1*x_n)+ c2 * r1* cos(r1*x_n);
    end
end


if(A2 == 0)
    du = du + ((bc2-bc1)/L);    
end

%particular part of the solution
if(k ~= 0)
    A = L^2/(A2*L^2-16*pi^2*k^2*A1);
    du = du+ A* 0.5 * k^2 * (4*pi*k/L) * cos((4*pi*k*x_n)/L);
end


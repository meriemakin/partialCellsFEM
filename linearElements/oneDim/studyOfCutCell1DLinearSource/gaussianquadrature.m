%---------------------------------------------------------------------------%                          
%                           Definition of the gauss points and weights      %
%                           for different accuracies                        %
%---------------------------------------------------------------------------%

%x_g and w are respectively the gauss point coordinate and the weight
%within a given accurate integration of a specific polynomial.
%j is the total number of gauss points and k is the gauss point number.
function [x_g, w] = gaussianquadrature(j,k);

%exact for polynomials of degree 3
if (j == 2)
    if(k==1)
         x_g = -sqrt(1/3);
         w = 1.0;
    elseif(k==2)
         x_g = sqrt(1/3);
         w = 1.0;
    end 
%exact for polynomials of degree 5
elseif (j == 3)
    if(k==1)
         x_g = -0.77459667;
         w = 0.55555555;
    elseif(k==2)
         x_g = 0.0;
         w = 0.88888889;
    elseif(k==3)
        x_g = 0.77459667;
         w = 0.55555555;
    end
%exact for polynomials of degree 7
elseif (j == 4)
    if(k==1)
         x_g = -0.86113631;
         w = 0.34785485;
    elseif(k==2)
         x_g = -0.33998104;
         w = 0.65214515;
    elseif(k==3)
        x_g = 0.33998104;
         w = 0.65214515;
    elseif(k==4)
        x_g = 0.86113631;
         w = 0.34785485;
    end    
end

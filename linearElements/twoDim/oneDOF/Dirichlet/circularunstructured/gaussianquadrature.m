%---------------------------------------------------------------------------%                          
%                           Definition of the gauss points and weights      %
%                           for different accuracies                        %
%---------------------------------------------------------------------------%

%xi_g and w1 are respectively the gauss point coordinate and the weight
%within a given accurate integration of a specific polynomial.
%eta_g and w2 are respectively the gauss point coordinate and the weight
%within a given accurate integration of a specific polynomial.
%j is the total number of gauss points and k is the gauss point number.
function [xi_g, eta_g, w1, w2] = gaussianquadrature(j,k,l);

%exact for polynomials of degree 3
if (j == 2)
    if(k==1)
         xi_g = -sqrt(1/3);
         w1 = 1.0;
    elseif(k==2)
         xi_g = sqrt(1/3);
         w1 = 1.0;
    end 
    if(l==1)
         eta_g = -sqrt(1/3);
         w2 = 1.0;
    elseif(l==2)
         eta_g = sqrt(1/3);
         w2 = 1.0;
    end 
%exact for polynomials of degree 5
elseif (j == 3)
    if(k==1)
         xi_g = -0.77459667;
         w1 = 0.55555555;
    elseif(k==2)
         xi_g = 0.0;
         w1 = 0.88888889;
    elseif(k==3)
        xi_g = 0.77459667;
         w1 = 0.55555555;
    end
    
    if(k==1)
         eta_g = -0.77459667;
         w2 = 0.55555555;
    elseif(k==2)
         eta_g = 0.0;
         w2 = 0.88888889;
    elseif(k==3)
        eta_g = 0.77459667;
         w2 = 0.55555555;
    end
%exact for polynomials of degree 7
elseif (j == 4)
    if(k==1)
         xi_g = -0.86113631;
         w1 = 0.34785485;
    elseif(k==2)
         xi_g = -0.33998104;
         w1 = 0.65214515;
    elseif(k==3)
        xi_g = 0.33998104;
         w1 = 0.65214515;
    elseif(k==4)
        xi_g = 0.86113631;
         w1 = 0.34785485;
    end    
    if(l==1)
         eta_g = -0.86113631;
         w2 = 0.34785485;
    elseif(l==2)
         eta_g = -0.33998104;
         w2 = 0.65214515;
    elseif(l==3)
        eta_g = 0.33998104;
         w2 = 0.65214515;
    elseif(l==4)
        eta_g = 0.86113631;
         w2 = 0.34785485;
    end    
%exact for polynomials of degree 9    
elseif (j == 5)
    if(k==1)
        xi_g = -0.90617985;
        w1 = 0.23692689;
    elseif(k==2)
        xi_g = -0.53846931;
        w1 = 0.47862867;
    elseif(k==3)
        xi_g = 0.0;
        w1 = 0.56888889;
    elseif(k==4)
        xi_g = 0.53846931;
        w1 = 0.47862867;
    elseif(k==5)
        xi_g = 0.90617985;
        w1 = 0.23692689;
    end    
    if(l==1)
        eta_g = -0.90617985;
        w2 = 0.23692689;
    elseif(l==2)
        eta_g = -0.53846931;
        w2 = 0.47862867;
    elseif(l==3)
        eta_g = 0.0;
        w2 = 0.56888889;
    elseif(l==4)
        eta_g = 0.53846931;
        w2 = 0.47862867;
    elseif(l==5)
        eta_g = 0.90617985;
        w2 = 0.23692689;
    end    
end

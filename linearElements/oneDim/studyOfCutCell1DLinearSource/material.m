%---------------------------------------------------------------------------%                          
%                           material property                               %
%---------------------------------------------------------------------------%

%x_n is the nodes coordinates
%i is the case number
function [A1, dA1] = material(xi, i);

%if the material property is the same for all elements
if (i == 0)
    A1 = 3.0;
    dA1 = 0.0;
%if the material property varies from element to element
elseif (i == 1)
end



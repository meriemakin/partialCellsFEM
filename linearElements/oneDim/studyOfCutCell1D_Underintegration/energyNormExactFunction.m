%---------------------------------------------------------------------------%                          
%                           Calculation of the energy norm                  %
%---------------------------------------------------------------------------%

%n is the number of gauss points
%elem is the element for which the norm is calculated
%K  the factor in the DE
%L the length of the system
%A1 factor in the DE
%A2 factor in the DE
%bc1 left BC
%bc2 right BC
function norm = energyNormExactFunction(n,elem, K,L,A1,A2, bc1, bc2);

norm = 0.0;
for i=1:n
  n
  [x_g, w] = gaussianquadrature(n,i);
  phi = shapefunction(x_g,length(elem));
  x = mapping(x_g,phi,elem);
  jacob = jacobian2(elem,x_g);
 %   x=0:0.0001:1;
    derivative = derivativeExactSolution(x, K,L,A1,A2, bc1, bc2);
 %   norm = trapz(x,derivative.*A1.*derivative);
  norm = norm + A1 * derivative^2 * w * (1/jacob);
end
    norm
    %norm = sqrt(norm);






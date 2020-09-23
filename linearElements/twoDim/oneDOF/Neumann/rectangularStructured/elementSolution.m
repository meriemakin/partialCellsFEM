%extracts the element nodal values from the global solution array using the
%element connectivity
function u_nE = elementSolution(solution, connEl);

[r,c] = size(connEl);

for i=1:c
    u_nE(i,1) = solution(connEl(i),1);
    
end

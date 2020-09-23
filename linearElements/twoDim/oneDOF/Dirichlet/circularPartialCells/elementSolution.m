%extracts the element nodal values from the global solution array using the
%element connectivity
function u_nE = elementSolution(solution, el);

[r,c] = size(el);
[r1,c1] = size(solution);

for i=1:c

    x = el(1,i);
    y = el(2,i);

    for j = 1:r1
        if(abs(x - solution(j,1)) < eps && abs(y - solution(j,2)) < eps)
             u_nE(i,1)  = solution(j,3);
             break;
        end
    end    
end

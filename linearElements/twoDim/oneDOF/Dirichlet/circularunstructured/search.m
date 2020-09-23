%---------------------------------------------------------------------------%                          
%                           generating nodes                                %
%---------------------------------------------------------------------------%

function [index] = search(edges, node1, node2);

[r,c] = size(edges);

for i=1:r
    if(edges(i,1) == node1 && edges(i,2)==node2)
       index  = i;
       return;
    elseif(edges(i,1) == node2 && edges(i,2) == node1)
       index  = i;
       return;
    end
end




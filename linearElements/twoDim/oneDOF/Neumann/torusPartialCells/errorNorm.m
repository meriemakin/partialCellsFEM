%---------------------------------------------------------------------------%                          
%                           Error norm                                      %
%---------------------------------------------------------------------------%

function [err,nin] = errorNorm(fem,analy);


[r,c] = size(analy);
[l1,l2] = size(fem);
err = 0;

for i=1:r
    x = analy(i,1)
    y = analy(i,2)
    sol = analy(i,3)
    
    for j=1:l1
       if(abs(fem(j,1) - x) <= 2*eps && abs(fem(j,2) - y) <= 2*eps)
            fem(j,1) 
            fem(j,2)
            err = err + abs(fem(j,3) - sol);
            break;
       end
    end
end


err = sqrt(err/r);
nin = r;



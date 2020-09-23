%---------------------------------------------------------------------------%                          
%                           Error norm                                      %
%---------------------------------------------------------------------------%

function [err,nin] = errorNorm(fem,analy);

nin = length(fem);

err = 0;

for j=1:nin
    err = err + abs(fem(j) - analy(j));
end

err = err/nin;



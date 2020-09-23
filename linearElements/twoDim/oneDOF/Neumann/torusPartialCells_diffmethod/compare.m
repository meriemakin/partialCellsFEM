%---------------------------------------------------------------------------%                          
%                           Compare the elements of two arrays              %
%---------------------------------------------------------------------------%

function truefalse = compare(indexes1, indexes2);

truefalse = 0;


for i=1:length(indexes1)
    el = indexes1(i);
    for j=1:length(indexes2)
        if( abs(el - indexes2(j)) <= 2*eps)
           truefalse = 1;
           return;
        end
    end
end

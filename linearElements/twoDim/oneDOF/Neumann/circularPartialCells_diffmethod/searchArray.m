%---------------------------------------------------------------------------%                          
%                           Searches Array for an Element                   % 
%---------------------------------------------------------------------------%

function index = searchArray(Array, el);

index = [];

if(~isempty(Array))
   for i=1:length(Array)
       if (abs(Array(i)- el) <= 2*eps)
           index = [index;i];
       end
   end
end

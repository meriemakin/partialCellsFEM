function [found] = searchindex(numelem,index)

found = 0;

for i=1:length(numelem)
	if(numelem(i,1) == index)
		found = 1;
		return;		
	end
end


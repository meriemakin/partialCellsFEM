%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%elements are the cartesian grid elements
%xC is the x coordinate of the center of the circular physical domain
%yC is the y coordinate of the center of the circular physical domain
%radius1 is the radius of the outer circular physical domain
%radius2 is the radius of the inner circular physical domain
function flags = elementflags(elems, xC, yC, radius1, radius2);

[n,r,c] = size(elems);
flags = zeros(c,1);


if (c == 1)
    flags = [2];
    return; 
end

%0 outside of domain
%1 totally inside of domain
%2 partially inside of domain
for i=1:c
    element = elems(:,:,i);
    count = 0;
    for j=1:r
    	x = element(1,j);
    	y = element(2,j);
        if((sqrt((xC- x)^2.0 + (yC - y)^2.0) < radius1) && (sqrt((xC- x)^2.0 + (yC - y)^2.0) > radius2) )
	  count = count + 1;
        end
    end
    if(count == 4)
      flags(i,1) = 1;
    elseif(count == 0)
      flags(i,1) = 0;
    elseif(count > 0 && count < 4)
      flags(i,1) = 2;
    end
end

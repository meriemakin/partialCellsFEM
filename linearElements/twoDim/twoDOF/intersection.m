%---------------------------------------------------------------------------%                          
%                           Definition of the mapping                       %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [xs,ys] = intersection(elem,xC, yC, radius);

xs = [];
ys = [];

x_n = elem(1,:);
y_n = elem(2,:);

minX = min(x_n);
maxX = max(x_n);
minY = min(y_n);
maxY = max(y_n);

%right element edge
y = yC + sqrt(radius^2 - (minX - xC)^2);
if(y >= minY && y <= maxY )
    if ((isequal(find(xs==minX), [])) && (isequal(find(ys==y),[])))
        xs = [xs;minX];
        ys = [ys;y];
    else
      if((find(xs == minX)) ~= (find(ys == y)))
    	xs = [xs;minX];
    	ys = [ys;y];
      end
    end
end 

%left element edge
y = yC + sqrt(radius^2 - (maxX - xC)^2);
if(y >= minY && y <= maxY )
    if ( (isequal(find(xs ==maxX), [])) && (isequal(find(ys==y), [])))     
        xs = [xs;maxX];
        ys = [ys;y];
    else
      if((find(xs == maxX)) ~= (find(ys == y)))
        xs = [xs;maxX];
        ys = [ys;y];
      end
    end

end

%lower element edge
x = xC + sqrt(radius^2 - (minY - yC)^2);
if(x >= minX && x <= maxX )
    if ( (isequal(find(xs ==x), [])) && (isequal((ys==minY), [])))     
        xs = [xs;x];
        ys = [ys;minY];
    else
      if((find(xs == x)) ~= (find(ys == minY)))
        xs = [xs;x];
        ys = [ys;minY];
      end
    end

end

%upper edge
x = xC + sqrt(radius^2 - (maxY - yC)^2);
if(x >= minX && x <= maxX )
    if ( (isequal(find(xs==x), [])) && (isequal(find(ys==maxY), [])))     
        xs = [xs;x];
        ys = [ys;maxY];
    else
      if((find(xs == x)) ~= (find(ys == maxY)))
        xs = [xs;x];
        ys = [ys;maxY];
      end
    end

end


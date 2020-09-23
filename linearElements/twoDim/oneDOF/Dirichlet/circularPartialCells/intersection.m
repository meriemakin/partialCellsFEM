%---------------------------------------------------------------------------%                          
%                           Intersection of physical domain with grid       %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%xC and yC are the coordinates of the center of the circular physical domain
%radius is the radius of the circular physical domain
function [xs,ys,edge,flags] = intersection(elem,xC, yC, radius);

xs = [];
ys = [];
flags = [];
edge = [];

x_n = elem(1,:);
y_n = elem(2,:);

minX = min(x_n);
maxX = max(x_n);
minY = min(y_n);
maxY = max(y_n);

%flags
% 0  intersection with no node
% 1  intersection with node


%edge
%1 left
%2 right
%3 lower
%4 upper

%left element edge

%fprintf('left element edge\n');
if(radius^2 - (minX - xC)^2 >= 0)
y1 = yC + sqrt(radius^2 - (minX - xC)^2);
y2 = yC - sqrt(radius^2 - (minX - xC)^2);

if(y1 >= (minY - 2*eps)  && y1 <= (maxY + 2*eps) )
    %if ((isequal(find(xs == minX), [])))
    %    if(y1 == minY || y1 == maxY)
    %       xs = [xs; minX];
    %       ys = [ys; y1];
    %       flags = [flags, 1];
    %    else
    %       xs = [xs; minX];
    %       ys = [ys; y1];
    %       flags = [flags; 0];
    %    end
    %else
      if((isempty(searchArray(xs,minX)) &&  isempty(searchArray(ys,y1))) || (compare(searchArray(xs,minX),searchArray(ys,y1)) == 0))
        if( abs(y1- minY) < 2*eps || abs(y1- maxY)<2*eps)
           xs = [xs; minX];
           ys = [ys; y1];
           flags = [flags; 1];
           edge = [edge; 1];
        else
           xs = [xs; minX];
           ys = [ys; y1];
           flags = [flags; 0];
           edge = [edge; 1];
        end
      end
    %end
elseif(y2 >= (minY - 2*eps)  && y2 <= (maxY + 2*eps) )
    %if ((isequal(find(xs==minX), [])))
    %  if(y2 == minY || y2 == maxY)
    %    xs = [xs; minX];
    %    ys = [ys; y2];
    %    flags = [flags; 1];
    %  else
    %    xs = [xs; minX];
    %    ys = [ys; y2];
    %    flags = [flags; 0];
    %  end
    %else
      
      if((isempty(searchArray(xs,minX)) &&  isempty(searchArray(ys,y2))) || (compare(searchArray(xs,minX),searchArray(ys,y2)) == 0)) 
       if( abs(y2- minY)<=2*eps || abs(y2- maxY)<=2*eps)
        xs = [xs; minX];
        ys = [ys;  y2];
        flags = [flags; 1];
        edge = [edge; 1];
       else
        xs = [xs; minX];
        ys = [ys; y2];
        flags = [flags; 0];
        edge = [edge; 1];
       end
      end
    %end

end
end
%xs
%ys 

%right element edge

%fprintf('right element edge\n');
if(radius^2 - (maxX - xC)^2 >= 0)
y1 = yC + sqrt(radius^2 - (maxX - xC)^2);
y2 = yC - sqrt(radius^2 - (maxX - xC)^2);
if(y1 >= (minY - 2*eps)  && y1 <= (maxY + 2*eps))
    %if ((isequal(find(xs ==maxX), []))) 
    %    if(y1 == minY || y1 == maxY)  
    %       xs = [xs; maxX];
    %       ys = [ys; y1];
    %       flags = [flags; 1];
    %    else
    %       xs = [xs; maxX];
    %       ys = [ys; y1];
    %       flags = [flags; 0];
    %    end
    %else
      if((isempty(searchArray(xs,maxX)) &&  isempty(searchArray(ys,y1))) || (compare(searchArray(xs,maxX),searchArray(ys,y1)) == 0)) 
        if(abs(y1- minY)<=2*eps || abs(y1- maxY)<=2*eps)
           xs = [xs; maxX];
           ys = [ys; y1];
           flags = [flags ; 1];
           edge = [edge; 2];
        else
           xs = [xs; maxX];
           ys = [ys; y1];
           flags = [flags; 0];
           edge = [edge; 2];
        end
      end
    %end
elseif(y2 >= (minY - 2*eps)  && y2 <= (maxY  + 2*eps) )
    %if ((isequal(find(xs == maxX), [])))
    %   fprintf('h2');
    %   if(y2 == minY || y2 == maxY)
    %      xs = [xs; maxX];
    %      ys = [ys; y2];
    %      flags = [flags; 1];
    %   else
    %      xs = [xs; maxX];
    %      ys = [ys; y2];
    %      flags = [flags; 0];
    %   end
    %else
      if((isempty(searchArray(xs,maxX)) &&  isempty(searchArray(ys,y2))) || (compare(searchArray(xs,maxX),searchArray(ys,y2)) == 0)) 
        if( abs(y2- minY)<=2*eps || abs(y2- maxY)<=2*eps)
           xs = [xs; maxX];
           ys = [ys; y2];
           flags = [flags; 1];
           edge = [edge; 2];
        else
           xs = [xs; maxX];
           ys = [ys; y2];
           flags = [flags; 0];
           edge = [edge; 2];
        end
      end
    %end
end
end
%xs
%ys

%lower element edge

%fprintf('lower element edge\n');

if(radius^2 - (minY - yC)^2 >= 0)
x1 = xC + sqrt(radius^2 - (minY - yC)^2);
x2 = xC - sqrt(radius^2 - (minY - yC)^2);
if(x1 >= (minX - 2*eps)  && x1 <= (maxX + 2*eps) )
    %if ((isequal(find(ys==minY), [])))     
    %    fprintf('hier');
    %    if(x1 == minX || x1 == maxX)
    %       %fprintf('hier');
    %       xs = [xs; x1];
    %       ys = [ys; minY];
    %       flags = [flags; 1];
    %    else
          %fprintf('hier');
    %       xs = [xs; x1];
    %       ys = [ys; minY];
    %       flags = [flags; 0];
    %    end
    %else
      if((isempty(searchArray(xs,x1)) &&  isempty(searchArray(ys,minY))) || (compare(searchArray(xs,x1),searchArray(ys, minY)) == 0))
        if(abs(x1- minX)<=2*eps || abs(x1-maxX)<=2*eps)
           xs = [xs; x1];
           ys = [ys; minY];
           flags = [flags; 1];
           edge = [edge; 3];
        else
           xs = [xs; x1];
           ys = [ys; minY]; 
           flags = [flags; 0];
           edge = [edge; 3];
        end
      end
     %end
elseif(x2 >= (minX - 2*eps) && x2 <= (maxX + 2*eps))
    %if ((isequal(find(ys==minY), []))) 
    %    if(x2 == minX || x2 == maxX)
    %       xs = [xs; x2];
    %       ys = [ys; minY];
    %       flags = [flags; 1];
    %    else
    %       xs = [xs; x2];
    %       ys = [ys; minY];
    %       flags = [flags; 0];
    %    end
    %else

      if((isempty(searchArray(xs,x2)) &&  isempty(searchArray(ys,minY))) || (compare(searchArray(xs,x2),searchArray(ys,minY)) == 0))
        if(abs(x2- minX)<=2*eps || abs(x2- maxX)<=2*eps)
           xs = [xs; x2];
           ys = [ys; minY];
           flags = [flags; 1];
           edge = [edge; 3];
        else
           xs = [xs; x2];
           ys = [ys; minY];
           flags = [flags; 0];
           edge = [edge; 3];
        end
      end
    %end

end
end
%xs
%ys

%upper edge

%fprintf('upper element edge\n');
if(radius^2 - (maxY - yC)^2 >= 0)
x1 = xC + sqrt(radius^2 - (maxY - yC)^2);
x2 = xC - sqrt(radius^2 - (maxY - yC)^2);

if(x1 >= (minX-2*eps)  && x1 <= (maxX+2*eps) )
    %if ((isequal(find(ys==maxY), [])))    
    %    if(x1 == minX || x1 == maxX)
    %       xs = [xs; x1];
    %       ys = [ys; maxY];
    %       flags = [flags; 1];
    %    else
    %       xs = [xs; x1];
    %       ys = [ys; maxY];
    %       flags = [flags; 0];
    %    end
    %else
      if((isempty(searchArray(xs,x1)) &&  isempty(searchArray(ys,maxY))) || (compare(searchArray(xs,x1),searchArray(ys,maxY)) == 0))
        if(abs(x1- minX)<=2*eps || abs(x1- maxX)<=2*eps)
           xs = [xs; x1];
           ys = [ys; maxY];
           flags = [flags; 1];
           edge = [edge; 4];
        else
           xs = [xs; x1];
           ys = [ys;  maxY];
           flags = [flags; 0];
           edge = [edge; 4];
        end
      end
    %end
elseif(x2 >= (minX - 2*eps)  && x2 <= (maxX + 2*eps) )
    %if ((isequal(find(ys==maxY), [])))
    %    if(x2 == minX || x2 == maxX)
    %       xs = [xs; x2];
    %       ys = [ys; maxY];
    %	   flags =  [flags ; 1];
    %    else
    %       xs = [xs; x2];
    %       ys = [ys; maxY];
    %       flags = [flags; 0];
    %    end
    %else
      if((isempty(searchArray(xs,x2)) &&  isempty(searchArray(ys,maxY))) || (compare(searchArray(xs,x2),searchArray(ys,maxY)) == 0))
        if( abs(x2-minX)<=2*eps || abs(x2- maxX)<=2*eps)
           xs = [xs; x2];
           ys = [ys; maxY];
           flags = [flags; 1];
           edge = [edge; 4];
        else
           xs = [xs; x2];
           ys = [ys; maxY];
           flags = [flags; 0];
           edge = [edge; 4];
        end
      end
    %end
end
end
%xs
%ys

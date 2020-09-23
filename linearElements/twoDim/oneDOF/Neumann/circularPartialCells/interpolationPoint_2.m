%---------------------------------------------------------------------------%                          
%                           Searches Array for an Element                   %
%---------------------------------------------------------------------------%

function [x,y,edge] = interpolationPoint_2(elmt, hnX, hnY, xs, ys, edge, flags, xC, yC, radius);


%edge
%1 left
%2 right
%3 lower
%4 upper

%find out edges where the hanging nodes belongs
if(hnX == elmt(1,1) && hnY == elmt(2,1))
   edge1 = 1;
   edge2 = 3;
elseif(hnX == elmt(1,2) && hnY == elmt(2,2))
   edge1 = 1;
   edge2 = 4;
elseif(hnX == elmt(1,3) && hnY == elmt(2,3))
   edge1 = 4;
   edge2 = 2;
elseif(hnX == elmt(1,4) && hnY == elmt(2,4))
   edge1 = 2;
   edge2 = 3;
end

%edge1
%edge2

%see if the hanging node is an intersection point itself
for i=1:length(xs)
    xt = xs(i);
    yt = ys(i);
    ed = edge(i);
    if( abs(xt - hnX) <= 2*eps && abs(yt - hnY)<= 2*eps && ed ~=0 )
%      fprintf('hanging node is an intersection point\n');
     	
        x_nE = elmt(1,:);
    	y_nE = elmt(2,:);

     	minX = min(x_nE);
     	maxX = max(x_nE);
     	minY = min(y_nE);
     	maxY = max(y_nE);
        
        if(abs(xt - minX) <=2*eps)
            xt = xt + ((maxX - minX)/5.0); 
	elseif(abs(xt-maxX) <=2*eps)
            xt = xt - ((maxX - minX)/5.0);
	end

        y1 = yC + sqrt(radius^2 - (xt - xC)^2);
        y2 = yC - sqrt(radius^2 - (xt - xC)^2);
        
        if(y1 > minY && y1 < maxY)
          yt = y1;
        elseif(y2 > minY && y2 < maxY)
          yt = y2;
        end
     	
        x = xt;
     	y = yt;
     	edge(i) = 0;
    	return;
    end
end

%see to which edge the intersection point belongs 
for i=1:length(xs)
    xt = xs(i);
    yt = ys(i);
    ed = edge(i);
    flag = flags(i);
    if((ed == edge1 || ed == edge2) && flag==0 && ed ~=0)
%       fprintf('second\n');
%       fprintf('which edge the intersection point');
       %check if this intersection point is not a node itself.
       
       %First Interpolation type: involves only one direction
       %x = xt;
       %y = yt;
       


       %Second interpolation: nodes involved
       %in the intepolation
       x = (xt + ((xs(1) + xs(2))*0.5))*0.5;
       
       y1 = yC + sqrt(radius^2 - (x - xC)^2);
       y2 = yC - sqrt(radius^2 - (x - xC)^2);
       miny = min(ys(1), ys(2));
       maxy = max(ys(1), ys(2));
       if(y1 > miny && y1 < maxy)
          y = y1;
       elseif(y2 > miny && y2 < maxy)
          y = y2;
       end

       edge(i) = 0;
       return;
    end    
end


%if I didnt find any collapsing intersection point or 
%an intersection point that belongs to a neighboring edge to
%the hanging node, create an interpolation point

%    fprintf('third\n');
%fprintf('no intersection no edge');
x =(xs(1) + xs(2))*0.5;
y1 = yC + sqrt(radius^2 - (x - xC)^2);
y2 = yC - sqrt(radius^2 - (x - xC)^2);
miny = min(ys(1), ys(2));
maxy = max(ys(1), ys(2));
if(y1 > miny && y1 < maxy)
   y = y1;
elseif(y2 > miny && y2 < maxy)
   y = y2;
end



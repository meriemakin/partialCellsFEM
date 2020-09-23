function len = curvelength(x1,y1,x2,y2,xC,yC,radius);

%find the corresonding quadrant

%first point
alpha = 0.0;

%first quad
if ( ((x1 - xC)>= 0.0) && ((y1 - yC)>=0) ) 
	alpha = acos((x1-xC)/radius);	
%second quad
elseif(((x1 -xC) < 0.0) && ((y1 - yC) > 0))
	alpha = acos((x1-xC)/radius);
%third quad
elseif(((x1 - xC) <= 0.0) && ((y1 -yC) <= 0.0))
	alpha = 2*pi - acos((x1-xC)/radius);
%fourth quad
elseif( ((x1 - xC) > 0.0) && ((y1 - yC) < 0.0))
	alpha = 2*pi - acos((x1-xC)/radius);
end

%second point
beta = 0.0;

%first quad
if ( ((x2 - xC)>= 0.0) && ((y2 - yC)>=0) ) 
        beta = acos((x2-xC)/radius);
%second quad
elseif(((x2 -xC) < 0.0) && ((y2 - yC) > 0))
        beta = acos((x2-xC)/radius);
%third quad
elseif(((x2 - xC) <= 0.0) && ((y2 -yC) <= 0.0))
        beta = 2*pi - acos((x2-xC)/radius);
%fourth quad
elseif( ((x2 - xC) > 0.0) && ((y2 - yC) < 0.0))
        beta = 2*pi - acos((x2-xC)/radius);
end

diff = abs(beta - alpha);

len = radius * diff;

function [xm,ym] = curvecenter(x1,y1,x2,y2,xC,yC,radius);


xm = 0.0;
ym = 0.0;


%midpoint of center of segment
xtil = 0.5*(x1 + x2);
ytil = 0.5*(y1 + y2);

%line passing through center of segment and center of circle
a = (yC - ytil)/(xC - xtil);
b = yC - a*xC;

%intersection of line with circle are generally two points
aa = 1 + a^2;
bb = -2*xC + 2*a*(b-yC);
cc = xC^2 + (b-yC)^2 - radius^2;


x1s = (-bb + sqrt(bb^2 - 4* aa*cc))/(2*aa);
x2s = (-bb - sqrt(bb^2 - 4* aa*cc))/(2*aa);

%corresponding y coordinates on circle

y1s = a*x1s + b;
y2s = a*x2s + b;


if((x1s >= min(x1,x2) && x1s <= max(x1,x2)) && (y1s >= min(y1,y2) && y1s <= max(y1,y2)))
	xm = x1s;
	ym = y1s;
elseif((x2s >= min(x1,x2) && x2s <= max(x1,x2)) && (y2s >= min(y1,y2) && y2s <= max(y1,y2)))
	xm = x2s;
	ym = y2s;
end

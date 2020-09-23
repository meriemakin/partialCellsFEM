function Dshape = shapeDerNeumann(x1,x2,y1,y2,x,y);


%   w2            w3
%   *-------------*
%   |             |
%   |             |
%   |             |
%   |             |
%   |             |
%   *-------------*
%   w1            w4

dw1dx = ((y2 - y)/(y2 - y1))* ((-1)/(x2 - x1));
dw1dy = ((-1)/(y2 - y1))* ((x2 - x)/(x2 - x1)); 
dw4dx = ((y2 - y)/(y2 - y1))* ((1)/(x2 - x1));
dw4dy = ((-1)/(y2 - y1))* ((x - x1)/(x2 - x1));
dw2dx = ((y - y1)/(y2 - y1))* ((-1)/(x2 - x1));
dw2dy = ((1)/(y2 - y1))* ((x2 - x)/(x2 - x1));
dw3dx = ((y - y1)/(y2 - y1))* ((1)/(x2 - x1));
dw3dy = ((1)/(y2 - y1))* ((x - x1)/(x2 - x1));

Dshape = [dw1dx dw2dx dw3dx dw4dx;
	  dw1dy dw2dy dw3dy dw4dy];



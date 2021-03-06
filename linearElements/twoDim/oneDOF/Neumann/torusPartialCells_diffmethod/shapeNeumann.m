%---------------------------------------------------------------------------%                          
%                           Definition of the shape function derivatives    %
%                           in the master element.                          %
%---------------------------------------------------------------------------%


function shape = shapeNeumann(x1,x2,y1,y2,x,y);

%   w2            w3
%   *-------------*
%   |             |
%   |             |
%   |             |
%   |             |
%   |             |
%   *-------------*
%   w1            w4
w1 = ((y2 - y)/(y2 - y1))* ((x2 - x)/(x2 - x1));
w4 = ((y2 - y)/(y2 - y1))* ((x - x1)/(x2 - x1));
w2 = ((y - y1)/(y2 - y1))* ((x2 - x)/(x2 - x1));
w3 = ((y - y1)/(y2 - y1))* ((x - x1)/(x2 - x1));

shape = [w1 w2 w3 w4];


function [a1, a2, a3, a4,rhs] = linearInterpolationRadiation(x1,x2,y1,y2,x,y,xC,yC,K,H,bc1,bc2);


nx = x - xC;
ny = y - yC;
len = sqrt(nx^2.0 + ny^2.0);
nx = nx / len;
ny = ny / len;
%len = sqrt(nx^2.0 + ny^2.0) 




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

dw1dx = ((y2 - y)/(y2 - y1))* ((-1)/(x2 - x1));
dw1dy = ((-1)/(y2 - y1))* ((x2 - x)/(x2 - x1));
dw4dx = ((y2 - y)/(y2 - y1))* ((1)/(x2 - x1));
dw4dy = ((-1)/(y2 - y1))* ((x - x1)/(x2 - x1));
dw2dx = ((y - y1)/(y2 - y1))* ((-1)/(x2 - x1));
dw2dy = ((1)/(y2 - y1))* ((x2 - x)/(x2 - x1));
dw3dx = ((y - y1)/(y2 - y1))* ((1)/(x2 - x1));
dw3dy = ((1)/(y2 - y1))* ((x - x1)/(x2 - x1));


a1 = (H*w1) + K*(dw1dx*nx + dw1dy * ny);
a2 = (H*w2) + K*(dw2dx*nx + dw2dy * ny);
a3 = (H*w3) + K*(dw3dx*nx + dw3dy * ny);
a4 = (H*w4) + K*(dw4dx*nx + dw4dy * ny);

rhs = bc2 + (H*bc1);



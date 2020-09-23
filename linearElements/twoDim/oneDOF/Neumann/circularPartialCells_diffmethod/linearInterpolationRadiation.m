function [w1, w2, w3, w4,rhs] = linearInterpolationRadiation(x1,x2,y1,y2,x,y,xC,yC,K,H,bc1,bc2);

x

y

nx = x - xC
ny = y - yC
len = sqrt(nx^2.0 + ny^2.0)
nx = nx / len
ny = ny / len
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

w1 = ((y2 - y)/(y2 - y1))* ((x2 - x)/(x2 - x1))
w4 = ((y2 - y)/(y2 - y1))* ((x - x1)/(x2 - x1))
w2 = ((y - y1)/(y2 - y1))* ((x2 - x)/(x2 - x1))
w3 = ((y - y1)/(y2 - y1))* ((x - x1)/(x2 - x1))

dw1dx = ((y2 - y)/(y2 - y1))* ((-1)/(x2 - x1))
dw1dy = ((-1)/(y2 - y1))* ((x2 - x)/(x2 - x1)) 
dw4dx = ((y2 - y)/(y2 - y1))* ((1)/(x2 - x1))
dw4dy = ((-1)/(y2 - y1))* ((x - x1)/(x2 - x1))
dw2dx = ((y - y1)/(y2 - y1))* ((-1)/(x2 - x1))
dw2dy = ((1)/(y2 - y1))* ((x2 - x)/(x2 - x1))
dw3dx = ((y - y1)/(y2 - y1))* ((1)/(x2 - x1))
dw3dy = ((1)/(y2 - y1))* ((x - x1)/(x2 - x1))


w1 = ((H/K)*w1) +  (dw1dx*nx + dw1dy * ny)
w2 = ((H/K)*w2) +  (dw2dx*nx + dw2dy * ny)
w3 = ((H/K)*w3) +  (dw3dx*nx + dw3dy * ny)
w4 = ((H/K)*w4) +  (dw4dx*nx + dw4dy * ny)

rhs = bc2 + ((H/K)*bc1)



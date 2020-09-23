function [normal] = normalNeumann(x,y,xC,yC);

nx = xC - x;
ny = yC - y;
len = sqrt(nx^2.0 + ny^2.0);
nx = nx / len;
ny = ny / len;

normal = [nx;ny];



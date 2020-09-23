%---------------------------------------------------------------------------%                          
%                           Definition of the mapping                       %
%---------------------------------------------------------------------------%

%x and y are the the real coordinates as a function of the isoparametric coordinate 
%phi are the shape functions
%element is the element in which x,y lies

function [x,y] = mapping(phi, elem);

xs = elem(1,:);
ys = elem(2,:);
x = phi*xs';
y = phi*ys';


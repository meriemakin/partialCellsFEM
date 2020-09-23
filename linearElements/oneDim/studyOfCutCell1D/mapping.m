%---------------------------------------------------------------------------%                          
%                           Definition of the mapping                       %
%---------------------------------------------------------------------------%

%x is the real coordinates as a function of the isoparametric coordinate 
%phi are the shape functions
%element is the element in which x,y lies

function [x] = mapping(phi, elem);

xs = elem(1,:);
x = phi*xs';


%---------------------------------------------------------------------------%                          
%                           Definition of the mapping                       %
%---------------------------------------------------------------------------%

%x is the real coordinate as a function of the isoparametric coordinate 
%xi is the isoparametric coordinate
%phi are the shape functions
%element is the element in which x lies
function x = mapping(xi,phi,elem);
x = 0;
for i=1:length(phi)
    x= x+ phi(i)*elem(1,i);
end

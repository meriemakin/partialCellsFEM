%---------------------------------------------------------------------------%                          
%                           Definition of the shape functio Matrix          %
%                           in the master element.                          %
%---------------------------------------------------------------------------%

%phi are the shape functions 
function matrix = shapefunctionMatrix(phi);

matrix = [phi(1) 0 phi(2) 0 phi(3) 0 phi(4) 0;
         0 phi(1) 0 phi(2) 0 phi(3) 0 phi(4)];

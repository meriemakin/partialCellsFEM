%---------------------------------------------------------------------------%                          
%                           Definition of the shape functio Matrix          %
%                           in the master element.                          %
%---------------------------------------------------------------------------%

%phi are the shape functions 
function matrix = BMatrix(DphiR);

matrix = [DphiR(1,1) 0 DphiR(2,1) 0 DphiR(3,1) 0 DphiR(4,1) 0;
         0 DphiR(1,2) 0 DphiR(2,2) 0 DphiR(3,2) 0 DphiR(4,2);
         DphiR(1,2) 0 DphiR(2,2) 0 DphiR(3,2) 0 DphiR(4,2) 0;
         0 DphiR(1,1) 0 DphiR(2,1) 0 DphiR(3,1) 0 DphiR(4,1)];

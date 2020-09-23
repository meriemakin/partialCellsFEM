%---------------------------------------------------------------------------%                          
%                           Definition of the shape function derivatives    %
%                           in the master element.                          %
%---------------------------------------------------------------------------%

%Dphi are the derivative of the shape functions with respect to the isoparametric coordinates
%F is the finite element deformation gradient (gradient of the actual coordinates with respect to the isoparametric coordinate system)
function DphiR = shapefunctionFirstDerivativeReal(Dphi, F);

Finv = F^-1;
xi_x = Finv(1,1);
xi_y = Finv(1,2);
eta_x = Finv(2,1);
eta_y = Finv(2,2);

DphiR(1,1) = Dphi(1,1) * xi_x + Dphi(1,2) * eta_x; 
DphiR(1,2) = Dphi(1,1) * xi_y + Dphi(1,2) * eta_y;
DphiR(2,1) = Dphi(2,1) * xi_x + Dphi(2,2) * eta_x;
DphiR(2,2) = Dphi(2,1) * xi_y + Dphi(2,2) * eta_y;
DphiR(3,1) = Dphi(3,1) * xi_x + Dphi(3,2) * eta_x;
DphiR(3,2) = Dphi(3,1) * xi_y + Dphi(3,2) * eta_y;
DphiR(4,1) = Dphi(4,1) * xi_x + Dphi(4,2) * eta_x;
DphiR(4,2) = Dphi(4,1) * xi_y + Dphi(4,2) * eta_y;


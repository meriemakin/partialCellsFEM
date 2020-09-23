%---------------------------------------------------------------------------%                          
%                           system matrix                                   %
%---------------------------------------------------------------------------%

function s = SystemMatrix(c_n,alpha, DeltaX);

nc = length(c_n);
deltaX = c_n(2,1) - c_n(1,1);
s = zeros(nc,nc);

for i=2:(nc-1)
    s(i,i-1) = 1;
    s(i,i) = -2;
    s(i,i+1) = 1; 
end


%first line ivolves the boundary condition
%and accounts for a ghost centroid passing the 
%left end of the domain

% ----x----[----x----|----x----|----x----|

s(1,1) = -4; 
s(1,2) = (4/3);


%first line involves the boundary condition
%and accounts for a ghost centroid passing the 
%left end of the domain

% ----x----|----x----]----x----

if(alpha == 1)
s(nc,nc-1)= (4/3);
s(nc,nc) = -4;
end


if(alpha < 1 && alpha ~=0.5)
x = (1 + alpha) * DeltaX;
N1 = ((x - (3*DeltaX/2))*(x - (5*DeltaX/2)))/(2*DeltaX^2)
N2 = -((x - (DeltaX/2))*(x - (5*DeltaX/2)))/(DeltaX^2)
N3 =  ((x - (DeltaX/2))*(x - (3*DeltaX/2)))/(2*DeltaX^2)
s(nc, nc-1) = (1 - (N1/N3));
s(nc,nc) = -(2 + (N2/N3));
end


%if we have a cut cell, but the boundary hits exactly the a centroid
%i.e. when alpha = 0.5
if( alpha == 0.5)
s(nc,nc) = 1;
end

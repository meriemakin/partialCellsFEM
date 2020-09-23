%---------------------------------------------------------------------------%                          
%                           system matrix                                   %
%---------------------------------------------------------------------------%

function s = SystemMatrix(c_n,alpha, DeltaX,flag);

nv = length(c_n)
s = zeros(nv,nv);

x = 0;
N2_x = (-2*x+1.5*DeltaX)/(0.5*DeltaX^2);
N3_x = (2*x-0.5*DeltaX)/(1.5*DeltaX^2);
s(1,1) = (-1/DeltaX)-N2_x;
s(1,2) = (1/DeltaX)-N3_x;


for i=2:nv-1
    s(i,i-1) = 1;
    s(i,i) = -2;
    s(i,i+1) = 1; 
end

if(alpha < 1)
if(flag == 0)
	x = (1.5+alpha)*DeltaX;
	N1_x = (2*x - (2.5 + alpha)*DeltaX)/((1.5+alpha)*DeltaX^2);
	N2_x = -(2*x - (1.5+alpha)*DeltaX)/((0.5+alpha)*DeltaX^2);
	s(nv, nv-2) = N1_x;
	s(nv, nv-1) = N2_x  + (1/DeltaX);
	s(nv, nv) = - (1/DeltaX);
elseif(flag == 1)
	s(nv,nv) = (1/DeltaX);
	s(nv,nv-1) = -(1/DeltaX);
end
end

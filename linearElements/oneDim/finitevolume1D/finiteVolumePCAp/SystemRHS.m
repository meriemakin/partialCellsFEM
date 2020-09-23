%---------------------------------------------------------------------------%                          
%                           Definition of the RHS of the system             %
%---------------------------------------------------------------------------%


function f = SystemRHS(DeltaX, alpha, K, A, nc, bc1, bc2,flag);


if(alpha >1)
  for i=1:20
    fprintf('The cut cell length is larger than the entire cell!!!\n');
  end
end

f = zeros(nc,1);

x = 0;
N1_x = (8*(x-DeltaX))/(3*DeltaX^2);
f(1,1) = ((-K/A)*DeltaX)+N1_x*bc1;

for i=2:nc-1
   f(i,1) = (-K*(DeltaX^2.0))/A;
end

%in the case of a cut cell
if(alpha < 1)
if(flag == 0)
	x = (1.5 + alpha) * DeltaX;
	N3_x = (2*x-DeltaX)/((alpha+1.5)*(alpha+0.5)*DeltaX^2);
	f(nc,1) =  (-(K/A)*alpha*DeltaX)- N3_x*bc2;
elseif(flag == 1)
	f(nc,1) = ((K/A)*alpha*DeltaX) + bc2;
end
end

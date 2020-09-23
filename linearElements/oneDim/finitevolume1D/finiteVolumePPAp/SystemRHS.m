%---------------------------------------------------------------------------%                          
%                           Definition of the RHS of the system             %
%---------------------------------------------------------------------------%


function f = SystemRHS(DeltaX, alpha, K, A, nc, bc1, bc2);


if(alpha >1)
  for i=1:20
    fprintf('The cut cell length is larger than the entire cell!!!\n');
  end
end
f = zeros(nc,1);
for i=1:nc
   f(i,1) = (-K*(DeltaX^2.0))/A;
end

%bc contributions
f(1,1)  =  f(1,1) + -(8/3)*bc1;

if(alpha == 1)
f(nc,1) =  f(nc,1) + (-8/3)*bc2;
end

%in the case of a cut cell
if(alpha < 1 && alpha ~=0.5)
%change only the contribution to f(nc,1)
x = (1 + alpha) * DeltaX;
N3 = ((x-(DeltaX/2))*(x-(3*DeltaX/2)))/(2*DeltaX^2);
f(nc,1) = f(nc,1) + (-1/N3)*bc2;
end


%if the cell end hits a centroid, i.e. alpha =0.5
if(alpha == 0.5)
f(nc,1) = bc2;
end

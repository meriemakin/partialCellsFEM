function flux = nodalForce(elem,solution,A,xC,yC, radius1, radius2);

% use the intersection point of the inner circular
%boundary with the element

[xs,ys,edge,flags] = intersection(elem,xC, yC,radius1,radius2);

x = xs(1,1); 
y = ys(1,1);

%find xi eta coordinates corresponding to a point
%luying on the inner circular boundary inside of the 
%considered element

x_nE = elem(1,:);
y_nE = elem(2,:);

minX = min(x_nE);
maxX = max(x_nE);
minY = min(y_nE);
maxY = max(y_nE);

aX = 2/(maxX - minX);
bX = 1 - aX * maxX;

aY = 2/(maxY - minY);
bY = 1 - aY * maxY;

xi = aX * x + bX;
eta = aY * y + bY;

%get the nodal solution to the element
sol = elementSolution(solution, elem)

%calculate the flux on that point
Dphi = shapefunctionFirstDerivative(xi,eta);
F = deformationGradient(elem,Dphi);
Finv = inv(F);

dudx = 0.0;
dudy = 0.0;

for i=1:4
    for j=1:2
    	dudx = dudx + Dphi(i,j)*Finv(j,1)*sol(i,1);
 	dudy = dudy + Dphi(i,j)*Finv(j,2)*sol(i,1);
    end
end

%calculate the normal direction to the inner circular
%boundary
nx = xC - x;
ny = yC - y;
len = sqrt(nx^2.0 + ny^2.0);
nx = nx / len;
ny = ny / len;

%get the flux value
flux = A * (dudx * nx + dudy * ny);


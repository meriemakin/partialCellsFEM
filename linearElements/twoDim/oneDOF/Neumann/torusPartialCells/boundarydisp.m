function disp = boundarydisp(elem,solution,A,xC,yC, radius1, radius2);

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
sol = elementSolution(solution, elem);

%nodal shape functions
phi = shapefunction(xi, eta);

disp = 0.0;

for i=1:4
    disp = disp + phi(i)*sol(i,1);
end


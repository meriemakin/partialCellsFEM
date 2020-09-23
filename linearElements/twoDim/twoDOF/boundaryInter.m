%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%n is the total number of gauss points used for integration
%nn is the number of nodes per element
%A is coefficient in differential equation to be solved
function [lhs,rhs] = boundaryInter(elem,bcflag,bc,xC,yC,radius);

lhs = [];
rhs = [];

if(bcflag == 1)
     [xs,ys] = intersection(elem,xC, yC, radius);
     x1 = min(elem(1,:));
     x2 = max(elem(1,:));
     y1 = min(elem(2,:));
     y2 = max(elem(2,:));
     lhs = zeros(length(xs)*2,8);
     rhs = zeros(length(xs)*2,1);
     for i=1:length(xs)
         [w1X, w2X, w3X, w4X,rhsX] = linearInterpolation(x1,x2,y1,y2,xs(i,1),ys(i,1),bc(1,1));
          lhs(2*(i-1) + 1,:)=[w1X 0 w2X 0 w3X 0 w4X 0];
          rhs(2*(i-1) + 1,1) = rhsX;
         [w1Y, w2Y, w3Y, w4Y,rhsY] = linearInterpolation(x1,x2,y1,y2,xs(i,1),ys(i,1),bc(2,1));
          lhs(2*(i-1) + 2,:)=[0 w1Y 0 w2Y 0 w3Y 0 w4Y];
          rhs(2*(i-1) + 2,1) = rhsY;
     end     
end

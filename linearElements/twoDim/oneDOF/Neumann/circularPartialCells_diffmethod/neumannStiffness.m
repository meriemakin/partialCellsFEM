%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%n is the total number of gauss points used for integration
%nn is the number of nodes per element
function s = neumannStiffness(x1,x2,y1,y2,nn,xs,ys,xC,yC,len);

s = zeros(nn,nn);

w = 0.5;

for i=1:2
	Dshape = shapeNeumann(x1,x2,y1,y2,xs(i),ys(i));
        normal = normalNeumann(xs(i),ys(i),xC,yC);	
	s = s + Dshape'*normal*normal'*Dshape*w*len;        
end


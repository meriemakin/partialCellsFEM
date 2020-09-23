%---------------------------------------------------------------------------%                          
%                           Definition of the element stiffness matrix      %
%---------------------------------------------------------------------------%

%elem is the element we are looking at
%n is the total number of gauss points used for integration
%nn is the number of nodes per element
function s = neumannStiffness(x1,x2,y1,y2,nn,xs,ys,xC,yC,len,radius,flag);

s = zeros(nn,nn);


if (flag == 1)
	w = 0.5;

	%trapezoidal rule
	for i=1:2
        	%shape = shapeNeumann(x1,x2,y1,y2,xs(i),ys(i));
		Dshape = shapeDerNeumann(x1,x2,y1,y2,xs(i),ys(i));
        	normal = normalNeumann(xs(i),ys(i),xC,yC);	
		%s = s + shape'*normal'*Dshape*w*len;  
		s = s + Dshape'*normal*normal'*Dshape*w*len;      
	end

elseif(flag == 2)
	%simpson's rule
	[xm,ym] = curvecenter(xs(1),ys(1),xs(2),ys(2),xC,yC,radius);
	xsn = [xs(1); xm; xs(2)];
	ysn = [ys(1); ym; ys(2)];


	weights  = [1/6; 4/6; 1/6];


	for i=1:3
        	Dshape = shapeDerNeumann(x1,x2,y1,y2,xsn(i),ysn(i));
        	normal = normalNeumann(xsn(i),ysn(i),xC,yC);
        	s = s + Dshape'*normal*normal'*Dshape*weights(i)*len;
	end

end

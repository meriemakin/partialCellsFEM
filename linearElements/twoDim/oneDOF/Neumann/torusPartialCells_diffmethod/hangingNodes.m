function hN = hangingNodes(elem,xC,yC,radius1,radius2,connElem);

hN = [];
[r,c] = size(elem);

for j=1:c
    x = elem(1,j);
    y = elem(2,j);
    dist = sqrt((xC- x)^2.0 + (yC - y)^2.0);
    if(((dist-radius1) >= 0.0) || ((dist-radius2) <= 0.0))
       if((abs(dist - radius1) - abs(dist - radius2)) < 0.0)
       	  hN=[hN;connElem(j),x, y,1];
       elseif((abs(dist - radius1)- abs(dist - radius2))> 0.0)
 	  hN=[hN;connElem(j),x, y,2];
       end
    end
end



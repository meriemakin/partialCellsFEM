%---------------------------------------------------------------------------%                          
%                           Assembling of the element stiffness matrix      %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

%elem are the elements
%conn is the connectivity
%ng is the number of gauss points in one direction
%rhs the value of the function f in the PDE
%tnn the total number of nodes in the domain
function [solution] = computeAndAssemble(elem,conn,ng,A,ff,tnn,flagsE,xC,yC,radius,bc1,bc2,H,x_n,y_n,Nx);

%allocate memory for the stiffness matrix
%create the stiffness matrix and the force vector only for the
%fully and partially covered elements. The elements that do not belong
%to the physical domain at all wont receive a contribution from the PDE 
[nelements, numberofnodes] = size(conn);
%s = zeros(tnn,tnn);
s = sparse(tnn, tnn);
f = zeros(tnn,1);

%make copies of the coordinates
x_nC = x_n;
y_nC = y_n;

for i = 1:nelements 

    if(flagsE(i) == 1 || flagsE(i) == 2 || flagsE(i) == 3)
    	sElem = elementStiffness(elem(:,:,i),ng,numberofnodes);
   
    	connElem = conn(i,:);
    
    	for j = 1:numberofnodes
        	for k = 1:numberofnodes
            		s(connElem(j),connElem(k)) = s(connElem(j),connElem(k)) + sElem(j,k);
        	end
    	end
    
    	fElem = elementForceVector(ng,numberofnodes,elem(:,:,i),A,ff);
    	for j = 1:numberofnodes
        	f(connElem(j),1) = f(connElem(j),1) + fElem(j,1);
    	end
     end

end


%delete the equations corresponding to the hanging nodes in the global
%system
for i=1:nelements

        connElem = conn(i,:);
        elemt = elem(:,:,i);

        if(flagsE(i) == 2 || flagsE(i) == 3)
                hN = hangingNodes(elemt,xC,yC,radius,connElem);
                [nr,nc] = size(hN);
                for j=1:nr
                        %can be made faster if flagging the nodes
                        %that have been seen already
                        s(hN(j,1),:)       = zeros(1,tnn);
                        f(hN(j,1))         = 0.0;
                end

        end


end

%
%The partially covered elements are subject to boundary conditions.
%Build the boundary conditions in the system
for i=1:nelements

     connElem = conn(i,:);
     elemt = elem(:,:,i);
    
     x_nE = elemt(1,:);
     y_nE = elemt(2,:);

     minX = min(x_nE);
     maxX = max(x_nE);
     minY = min(y_nE);
     maxY = max(y_nE);
     
     %if element partially covered by physical domain
     if(flagsE(i) == 2)
        
        %fprintf('element %6.2f \n',i);       
 
     	%get the hanging nodes
     	hN = hangingNodes(elemt,xC,yC,radius,connElem); 
     	[nr,nc] = size(hN);
     
     	%get the intersection coordinates of the physical domain with 
     	%the rectangular element
     	[xs,ys,edge,flags] = intersection(elemt,xC, yC, radius);
     	ni = length(xs);
 
     	for j=1:nr
        	%get the corresponding intersection point to the hanging node
         	[x,y,edge] = interpolationPoint(elemt, hN(j,2), hN(j,3), xs, ys, edge, flags, xC, yC, radius);

		[w1, w2, w3, w4,rhs] = linearInterpolationRadiation(minX,maxX,minY,maxY,x,y,xC,yC,A,H,bc1,bc2);

               	s(hN(j,1),connElem(1)) = s(hN(j,1),connElem(1)) + w1;
               	s(hN(j,1),connElem(2)) = s(hN(j,1),connElem(2)) + w2;
               	s(hN(j,1),connElem(3)) = s(hN(j,1),connElem(3)) + w3;
               	s(hN(j,1),connElem(4)) = s(hN(j,1),connElem(4)) + w4;
               	f(hN(j,1)) = f(hN(j,1)) + rhs;
      	end
      end

end

%Delete all entries in the system that are associated with nodes
%that do belong to the outside of domain elements

iterate = 1;
pos = 1;

while (iterate == 1) 

     %check if the diagonal entry is different from zero
     %if it is equal to zero, then it must be the line associted with
     %a node that does belong to an outside of domain element 
     if(abs(s(pos,pos) - 0.0) <= 2*eps)
        %So delete rows and columns associated with them
        s(pos,:) = [];
        f(pos) = [];
        s(:,pos) = [];
        x_nC(pos) = [];
        y_nC(pos) = [];
     else
         pos = pos + 1;
     end

     if(pos > length(f))
        iterate = 0;
        break;
     end

end


%solve system of equations
u_n = s\f;


solution = [x_nC, y_nC, u_n];


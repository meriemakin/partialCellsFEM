%---------------------------------------------------------------------------%                          
%                           Assembling of the element stiffness matrix      %
%                           and the element force vector                    %
%---------------------------------------------------------------------------%

%elem are the elements
%conn is the connectivity
%ng is the number of gauss points in one direction
%rhs the value of the function f in the PDE
%tnn the total number of nodes in the domain
function [solution] = computeAndAssemble(elem,conn,ng,A,ff,tnn,flagsE,xC,yC,radius,bc1,bc2,x_n,y_n,Nx);

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

%(s-s')
%s
%f

epsilon = 1000;
dof     = ((Nx + 1)*(Nx/2)) + (Nx/2) + 1;
s(dof, dof) = s(dof,dof)  + epsilon;
f(dof)  = f(dof) + epsilon * bc1; 

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


%f

%dof = ((Nx + 1)*(Nx/2)) + (Nx/2) + 1;
%middlex = x_nC(dof);
%middley = y_nC(dof);
%[s,f] = applyBCDirichlet(s,f,bc1,dof);
%x_nC(dof) = [];
%y_nC(dof) = [];


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
%               hN(j,2)
%               hN(j,3)
        	%get the corresponding intersection point to the hanging node
         	[x,y,edge] = interpolationPoint(elemt, hN(j,2), hN(j,3), xs, ys, edge, flags, xC, yC, radius);

          	%write an interpolation equation per hanging node
          	%s(hN(j,1),:)       = zeros(1,tnn);
                %s
                %f
          	%[w1, w2, w3, w4,rhs] = linearInterpolation(minX,maxX,minY,maxY,x,y,bc1);
		[w1, w2, w3, w4,rhs] = linearInterpolationNeumann(minX,maxX,minY,maxY,x,y,xC,yC,A,bc2);

               	s(hN(j,1),connElem(1)) = s(hN(j,1),connElem(1)) + w1;
               	s(hN(j,1),connElem(2)) = s(hN(j,1),connElem(2)) + w2;
               	s(hN(j,1),connElem(3)) = s(hN(j,1),connElem(3)) + w3;
               	s(hN(j,1),connElem(4)) = s(hN(j,1),connElem(4)) + w4;
               	f(hN(j,1)) = f(hN(j,1)) + rhs;
		%hN 
                %f
                %pause
      	end
      end

end

%Delete all entries in the system that are associated with nodes
%that do belong to the outside of domain elements

%s
%f


%dof = ((Nx + 1)*(Nx/2)) + (Nx/2) + 1;

%get elements to which the middle node belongs
%middleel = [ Nx * ((Nx * 0.5) -1) + (Nx * 0.5), Nx * ((Nx * 0.5) -1) + (Nx * 0.5) + 1, Nx * (Nx * 0.5) + (Nx * 0.5), Nx * (Nx * 0.5) + (Nx * 0.5) + 1];

%s(dof,:)       = zeros(1,tnn);
%f(dof)         = 0.0;

%s(dof,dof) = s(dof,dof) +  1.0;
%f(dof) = f(dof) + bc1;

%s
%f

%for i=1:4
%      me = middleel(i);
%end
%middlex = x_nC(dof);
%middley = y_nC(dof);

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

%s

%enforce a dirichlet boundary condition at the center of the circular 
%domain. The mesh is so set up such that a node falls always on the
%center of the circle

%dof = 0;

%for count=1:length(x_nC)

%	if(abs(x_nC(count)- middlex)<eps && abs(y_nC(count)- middley)<eps)
%		dof = count;
%	end
%end

%dof

%degree of freedom of the center node
%dof = ((Nx + 1)*(Nx/2)) + (Nx/2) + 1;
%[s,f] = applyBCDirichlet(s,f,bc1,dof);
%s
%f

%s

%x_nC(dof) = [];
%y_nC(dof) = [];

%solve system of equations
u_n = s\f;

%insert dirichlet bc
%u_n = [u_n;bc1];
%x_nC = [x_nC;middlex];
%y_nC = [y_nC;middley];

%for i=1:length(x_nC)
%    if(s(i,i) == 0)
%        i
%    end
%end

solution = [x_nC, y_nC, u_n];


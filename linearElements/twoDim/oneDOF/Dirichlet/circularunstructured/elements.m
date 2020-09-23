%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%xC and yC are the coordinates of the circle
%R is the radius of the circle
%rl is the refinement level
%elem are the current elements of the discretized domain
%conn is the connectivity of all elements
%h is the element size along the horizontal diameter
%flags stores the global dof of the boundary nodes
function [elem,conn,h,flags,nn,x_n,y_n,nel] = elements(xC,yC,R,rl);


%master nodes
theta = pi/4;

x_n = [xC - R, xC - R*cos(theta), xC, xC + R*cos(theta), xC + R, xC + R*cos(theta), xC, xC - R*cos(theta), xC - (0.5*R*cos(theta)), xC, xC + (0.5*R*cos(theta)), xC+0.5*R, xC+(0.5*R*cos(theta)), xC, xC - (0.5*R*cos(theta)), xC - 0.5*R, xC];

y_n = [yC, yC - R*sin(theta), yC - R, yC - R*sin(theta), yC, yC + R*sin(theta), yC + R, yC + R*sin(theta), yC - (0.5*R*sin(theta)), yC - 0.5*R, yC - 0.5*R*sin(theta), yC, yC + 0.5*R*sin(theta), yC + 0.5*R, yC + 0.5*R*sin(theta), yC, yC ];

conn = [1 2 9 16; 
        2 3 10 9;
        3 4 11 10;
        4 5 12 11;
        5 6 13 12;
        6 7 14 13;
        7 8 15 14;
        8 1 16 15;
        9 10 17 16;
        10 11 12 17;
        17 12 13 14;
        16 17 14 15];

%flag the nodes that lie on the outer boundary of the domain
flags = [1 2 3 4 5 6 7 8];

%edges of the entire mesh
edges = [1 2; 2 9; 9 16; 16 1; 2 3 ; 3 10; 10 9; 3 4 ; 4 11; 11 10; 4 5; 5 12; 12 11; 5 6; 6 13; 13 12; 6 7; 7 14;7 8; 14 13; 8 15; 15 14; 8 1; 15 16;16 17; 17 12; 17 14; 10 17];
%an array for storing the new numbering of the midpoints of the edges
midpoints = sparse(1,28);

%create the elements out of these nodes
for i =1:12
     elem(:,:,i) = [x_n(conn(i,1)),x_n(conn(i,2)),x_n(conn(i,3)),x_n(conn(i,4));
                   y_n(conn(i,1)),y_n(conn(i,2)),y_n(conn(i,3)),y_n(conn(i,4))];
end
%number of nodes and elements for the coarsest mesh
nn = length(x_n);
nel = 12;

%if one requested the coarsest mesh, then we are done
if(rl == 1)
  for i=1:nel
    for j=1:4
        plot(elem(1,j,i),elem(2,j,i),'+');
        hold on
    end
  end
   h = abs(elem(1,1,1) - elem(1,4,1));
   return;
%otherwise the coarsest has to be refined until the requested refinement level
else
   for k=2:rl
       
       %storage for the finer mesh
       elemc = [];
       connc = [];
       edgesc = [];
       nelc = 0;
       nnc = nn;

       %go over all elements in the coarse mesh
       for j=1:nel
              
         %take each element in the coarser mesh
         elemt = elem(:,:,j);
         
         %The way I implemented the data structure is that if the element has an edge on the boundary, the first two nodes in the data structure lie on the boundary. So we only need to check if the first node has a flag to be on the boundary 

	%splitting of non-boundary edges happen in the regular way
        x_nm = [];
        y_nm = []; 
        dofs = [];
        if(isempty(find(flags == conn(j,1))))        
         %create midpoints on the edges and a centroid
         x_nm = [0.5*(elemt(1,1)+elemt(1,2)),0.5*(elemt(1,2)+elemt(1,3)),0.5*(elemt(1,3)+elemt(1,4)),0.5*(elemt(1,4)+elemt(1,1)),0.25*(elemt(1,1)+elemt(1,2)+elemt(1,3)+elemt(1,4))]; 
         y_nm = [0.5*(elemt(2,1)+elemt(2,2)),0.5*(elemt(2,2)+elemt(2,3)),0.5*(elemt(2,3)+elemt(2,4)),0.5*(elemt(2,4)+elemt(2,1)),0.25*(elemt(2,1)+elemt(2,2)+elemt(2,3)+elemt(2,4))]; 
         %splitting of the edges that do lie on the boundary happen
         % in a different manner because one wants to keep the 
         % shape of a circle as most as possible, i.e. one wants the middle
         % point of the boundary edge to lie on the circle boundary!
         else
         xm  = 0.5*(elemt(1,1) + elemt(1,2));
         ym  = 0.5*(elemt(2,1) + elemt(2,2));
         sl = (ym - yC)/(xm - xC);
         int = ym - sl * xm;
         a = (1 + sl^2);
         b = (-2*xC) + 2 * sl * (int - yC);
         c = (-R^2) + xC^2 + (int-yC)^2;
         x1 = (-b  + sqrt(b^2 - 4 * a * c))/(2*a);
         x2 = (-b  - sqrt(b^2 - 4 * a * c))/(2*a);
         %dx = xC - xm;
         %dy = yC - ym;
         %dr = sqrt(dx^2 + dy^2);
         %dd = xm * yC - xC * ym;
          
         x_nm1  = [];
         if(((x1-elemt(1,2))<=0.0 && (x1-elemt(1,1))>=0.0) || ((x1-elemt(1,1))<=0.0 && (x1-elemt(1,2))>=0.0))
            x_nm1 = x1;     
         elseif(((x2-elemt(1,1))>=0.0 && (x2-elemt(1,2))<=0.0) || ((x2-elemt(1,2))>=0.0 && (x2-elemt(1,1))<=0.0))
            x_nm1 = x2;
         end
         x_nm = [x_nm1,0.5*(elemt(1,2)+elemt(1,3)),0.5*(elemt(1,3)+elemt(1,4)),0.5*(elemt(1,4)+elemt(1,1)),0.25*(elemt(1,1)+elemt(1,2)+elemt(1,3)+elemt(1,4))]; 
         y_nm = [sl*x_nm1+int,0.5*(elemt(2,2)+elemt(2,3)),0.5*(elemt(2,3)+elemt(2,4)),0.5*(elemt(2,4)+elemt(2,1)),0.25*(elemt(2,1)+elemt(2,2)+elemt(2,3)+elemt(2,4))]; 
         end

         %create 4 elements out of each element and save its connectivity
	 %edge one
         [index] = search(edges, conn(j,1), conn(j,2));

         if((midpoints(index) == 0))
            nnc = nnc + 1;
            dofs = [dofs, nnc];
            midpoints(index) = nnc;
            x_n = [x_n, x_nm(1)];
            y_n = [y_n, y_nm(1)];
         else
            dofs = [midpoints(index)];
         end

         [index ] = search(edges, conn(j,2), conn(j,3));

         if((midpoints(1,index) == 0))
            nnc = nnc + 1;
            dofs = [dofs, nnc];
            midpoints(1,index) = nnc;
            x_n = [x_n, x_nm(2)];
            y_n = [y_n, y_nm(2)];
         else
            dofs = [dofs, midpoints(1,index)];
         end
         
         [index ] = search(edges, conn(j,3), conn(j,4));

         if((midpoints(1,index) == 0))
            nnc = nnc + 1;
            dofs = [dofs, nnc];
            midpoints(1,index) = nnc;
            x_n = [x_n, x_nm(3)];
            y_n = [y_n, y_nm(3)];
         else
            dofs = [dofs, midpoints(1,index)];
         end

         [index ] = search(edges, conn(j,4), conn(j,1));

         if((midpoints(1,index) == 0))
            nnc = nnc + 1;
            dofs = [dofs, nnc];
            midpoints(index) = nnc;
            x_n = [x_n, x_nm(4)];
            y_n = [y_n, y_nm(4)];
         else
            dofs = [dofs, midpoints(1,index)];
         end

         nnc = nnc + 1;
         x_n = [x_n, x_nm(5)];
         y_n = [y_n, y_nm(5)];
 
         if(~(isempty(find(flags == conn(j,1)))))
            flags = [flags, dofs(1)];
         end 
        
         elemc(:,:,end+1) = [elemt(1,1),x_nm(1),x_nm(5),x_nm(4);
                            elemt(2,1),y_nm(1),y_nm(5),y_nm(4)];
         if(j == 1)
            elemc(:,:,1) = [];
         end

         connc(end+1,:) = [conn(j,1),dofs(1),nnc,dofs(4)]; 
         edgesc = [edgesc; conn(j,1),dofs(1); dofs(1),nnc; nnc,dofs(4); dofs(4), conn(j,1)];

         
         elemc(:,:,end+1) = [x_nm(1),elemt(1,2),x_nm(2),x_nm(5);
                            y_nm(1),elemt(2,2),y_nm(2),y_nm(5)];
         
         connc(end+1,:) = [dofs(1),conn(j,2),dofs(2),nnc];
         edgesc = [edgesc; dofs(1),conn(j,2); conn(j,2),dofs(2); dofs(2),nnc; nnc, dofs(1)];

         elemc(:,:,end+1) = [x_nm(5),x_nm(2),elemt(1,3),x_nm(3);
                             y_nm(5),y_nm(2),elemt(2,3),y_nm(3)];

         
         connc(end+1,:) = [nnc,dofs(2),conn(j,3),dofs(3)];
         edgesc = [edgesc; nnc,dofs(2); dofs(2),conn(j,3); conn(j,3),dofs(3); dofs(3), nnc];

         elemc(:,:,end+1) = [x_nm(4),x_nm(5),x_nm(3),elemt(1,4);
                             y_nm(4),y_nm(5),y_nm(3),elemt(2,4)];

         
         connc(end+1,:) = [dofs(4),nnc,dofs(3),conn(j,4)];
         edgesc = [edgesc; dofs(4),nnc; nnc,dofs(3); dofs(3),conn(j,4); conn(j,4), dofs(4)];
         
         % we substitute each element by 4 elements
         nelc = nelc + 4;

      end
 
         %update the mesh   
         nn = nnc;
         nel = nelc;
         elem = elemc;
         conn = connc;
         edges = edgesc;
         midpoints = sparse(1,length(edgesc));
   end
end


h = abs(elem(1,1,1) - elem(1,4,1));

%figure
%plot(x_n, y_n, '+')

%figure
%for i=1:nel
%    %elem(1,:,i)
%    %elem(2,:,i)
%    for j=1:4
%        plot(elem(1,j,i),elem(2,j,i),'o');
%        hold on
%    end
%end


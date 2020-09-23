%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history

tic

clc
clear all

%Length
Lx = 2;

%width
Ly = Lx;

%Number of elements in x-direction
Nx = 32; 
%Number of elements in y-direction
Ny = Nx;


%number of gauss points in one direction
ng = 2;

%coordinate of center of circular physical domain
xC = Lx / 2.0;
yC = Ly / 2.0;

%radius of outer circular physical domain
radius1 = (Lx/2.0);

%raidus of inner circle
radius2 = (Lx/4.0);

%Material factor in the DE
A = 1.0;

%The rhs of the system
rhs = 0.0;

%BCs on the boundary of the circular physical domain 

%dirichlet
bc1 = 0.5;
%dirichlet
bc2 = 0.862942655531861;


%create nodes
[x_n,y_n] = nodes(Lx, Ly, Nx, Ny);
plot(x_n, y_n, 'o');
hold on;
Center = [xC,yC];
circle(Center,radius1,1000,'b-');
circle(Center,radius2,1000,'b-');
drawnow;
%pause

%total number of nodes
tnn = length(x_n);

%the connectivity
conn = connectivity(Nx, Ny,2,2);

%create elements
elem = elements(x_n, y_n,conn);

%flag elements whether they are inside of physical domain, outside of physical domain or partially covered by physical domain.
flagsE = elementflags(elem, xC, yC, radius1,radius2);


% compute the element stiffness matrix and force vector
%and assemble them
[solution] = computeAndAssemble(elem, conn,ng,A,rhs,tnn,flagsE,xC,yC,radius1,radius2,bc1,bc2,x_n,y_n);
figure
plot3(solution(:,1),solution(:,2),solution(:,3),'+');
%calculate the analytical solution
%analytical solution does not involve a source term 
[y] = analyticalsolution(x_n,y_n,xC,yC,radius1,radius2,bc1,bc2);
hold on
plot3(y(:,1),y(:,2),y(:,3),'or')
%errorNorm
%[err, nin] = errorNorm(solution,y);

%flag = 2; 
%err =  H1errornorm(ng,conn,elem,solution,flag,A,bc1,bc2,xC,yC,radius1,radius2,flagsE)
%h = x_n(2)-x_n(1)

flux = nodalForce(elem(:,:,270),solution,A,xC,yC, radius1, radius2);
%some results

%H0-error
%h = [0.0078125, 0.015625000000000,0.031250000000000, 0.062500000000000,0.125000000000000,0.250000000000000 ];
%h0err = [3.456431878381728e-06,1.744737910807809e-05,7.395721797765163e-05, 3.926320167424314e-04, 0.001286334081929, 0.006315631269809];
%h1err= [0.005018296495402,0.010495592440974,0.021538325761759,0.043673060913271,0.083636616153183,0.176029290992293];
%enerr=[0.005018295305066,0.010495577939113,0.021538198786100,0.043671295946350,0.083626723636553,0.175915957459559];

toc

%data1 = solution;
%x1=data1(:,1);
%y1=data1(:,2);
%z1=data1(:,3);
%a1=size(data1);
%b1=a1(:,1);
%xlin1=linspace(min(x1),max(x1),b1);
%ylin1=linspace(min(y1),max(y1),b1);
%[X1,Y1]=meshgrid(xlin1,ylin1);
%Z1=griddata(x1,y1,z1,X1,Y1);

%data2 = y;
%x2=data2(:,1);
%y2=data2(:,2);
%z2=data2(:,3);
%a2=size(data2);
%b2=a2(:,1);
%xlin2=linspace(min(x2),max(x2),b2);
%ylin2=linspace(min(y2),max(y2),b2);
%[X2,Y2]=meshgrid(xlin2,ylin2);
%Z2=griddata(x2,y2,z2,X2,Y2);

%figure
%contourf(X1,Y1,Z1);
%figure
%contourf(X2,Y2,Z2)

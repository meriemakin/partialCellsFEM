%clear command window and variable history

tic

clc
clear all

format long

%Length
Lx = 2;

%width
Ly = Lx;

%Number of elements in x-direction
Nx = 64;

%Number of elements in y-direction
Ny = Nx;

%number of gauss points in one direction
ng = 2;

%coordinate of center of circular physical domain
xC = Lx / 2.0;
yC = Ly / 2.0;

%radius of circular physical domain
radius = (Lx/2.0);

%source term
rhs  = 2.0;

%material term (conductivity)
A = 1.0;

%Dirichlet BC (surface and surrounding medium temperature) 
bc1 = 2.0;

%Flux
bc2 = 1.0;

%Radiation
H = 0.5;

%create nodes
[x_n,y_n] = nodes(Lx, Ly, Nx, Ny);
plot(x_n, y_n, 'o');
hold on;
Center = [xC,yC];
circle(Center,radius,1000,'b-');
drawnow;
%pause

%total number of nodes
tnn = length(x_n);

%the connectivity
conn = connectivity(Nx, Ny,2,2);

%create elements
elem = elements(x_n, y_n,conn);

%flag elements whether they are inside of physical domain, outside of physical domain or partially covered by physical domain.
flagsE = elementflags(elem, xC, yC, radius);

% compute the element stiffness matrix and force vector
%and assemble them
[solution] = computeAndAssemble(elem,conn,ng,A,rhs,tnn,flagsE,xC,yC,radius,bc1,bc2,H,x_n,y_n,Nx);

%calculate the analytical solution
[y] = analyticalsolution(rhs,A,x_n,y_n,xC,yC,radius,bc1,bc2,H);
figure
plot3(y(:,1),y(:,2),y(:,3),'ro');
hold on
plot3(solution(:,1),solution(:,2),solution(:,3),'b+');
drawnow
%errorNorm
%[err, nin] = errorNorm(solution,y);

flag = 0; 
err =  H1errornorm(ng,conn,elem,solution,flag,rhs,A,bc1,bc2,xC,yC,radius,flagsE,H)
h = x_n(2)-x_n(1)

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
                     
%Results


%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
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
Nx = 8; 
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
bc1 = 3.0;

%neumann
bc2 = 0.6;

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
[y] = analyticalsolution(A,x_n,y_n,xC,yC,radius1,radius2,bc1,bc2,rhs);
hold on;
plot3(y(:,1),y(:,2),y(:,3),'or');
%errorNorm
%[err, nin] = errorNorm(solution,y);

% 8 20
% 16 70
%%32 270
%disp = boundarydisp(elem(:,:,270),solution,A,xC,yC, radius1, radius2);

flag = 1; 
err =  H1errornorm(ng,conn,elem,solution,flag,A,bc1,bc2,xC,yC,radius1,radius2,flagsE,rhs);
h = x_n(2)-x_n(1);

%some results

%H0-error
%h = [0.25,0.125,0.0625,0.03125,0.015625,0.0078125,0.00390625];
%h0error = [0.006489465886404,0.002058427606645,4.454760168961685e-04,1.670113136217037e-04,1.190168086179875e-04,4.034345852147729e-05,1.817685630424053e-05];
%h1error = [0.078452845027203,0.038482660108002,0.017837167687853,0.008718701612569,0.004366385449263,0.002131536312066,0.001061788164841];
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


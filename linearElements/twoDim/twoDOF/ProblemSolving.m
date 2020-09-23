%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history
clc
clear all


%Rectangular Geometry

%Length
Lx = 2;

%Width
Ly = Lx;

%number of nodes per element
nn = 4;

%Number of elements in x-direction
Nx = 2;

%Number of elements in y-direction
Ny = 2;

%penalty number
pen = 1000;

%temperature bc
T0 = 150;

%number of gauss points
ng = 4;

%coordinate of center of circular physical domain
xC = Lx/2.0;
yC = Ly/2.0;

%radius of circular physical domain
radius = (Lx/2.0); 

%Material factor in the DE
A = 1.0;

%The rhs of the system
rhs = [1.0; 1.0];

%BCs on upper right quarter of the circular physical domain 
bc = [0; 0];

%create nodes
[x_n,y_n] = nodes(Lx, Ly, Nx, Ny);
plot(x_n, y_n, 'o');
hold on;
Center = [xC,yC];
circle(Center,radius,1000,'b-');

%the connectivity
conn = connectivity(Nx, Ny,2,2);

%create elements
elements = elements(x_n, y_n,conn);

%flag elements whether they are inside of physical domain, outside of physical domain or partially covered by physical domain.
flagsE = elementflags(elements, xC, yC, radius);

%degree of freedom connectivity
dofconn = dofconnectivity(conn);

%bcflags
bcflags = bcflags(elements,flagsE,xC,yC, radius);

% compute the element stiffness matrix and force vector
%and assemble them
[s,f] = computeAndAssemble(elements,flagsE,bcflags, dofconn,length(x_n),ng,rhs,A,bc,xC,yC,radius);

%solve problem
u_n = s\f;

newx_n = zeros(length(x_n),1);
newy_n = zeros(length(y_n),1);

%plot solution
for i=1:length(x_n)
	newx_n(i,1) = x_n(i,1) + u_n(2*(i-1) + 1 ,1);
	newy_n(i,1) = y_n(i,1) + u_n(2*(i-1) + 2, 1);
end

hold on
plot(newx_n, newy_n, 'r*');


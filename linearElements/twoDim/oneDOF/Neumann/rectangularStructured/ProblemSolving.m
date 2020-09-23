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
Nx = 2; 

%Number of elements in y-direction
Ny = Nx;


%number of gauss points in one direction
ng = 2;


%Material factor in the DE
A = 1.0;

%The rhs of the system
rhs = 0.0;

%BCs  
bc = 2.0;

%create nodes
[x_n,y_n] = nodes(Lx, Ly, Nx, Ny);
%plot(x_n, y_n, 'o');

%total number of nodes
tnn = length(x_n);

%the connectivity
conn = connectivity(Nx, Ny,2,2);

%create elements
elem = elements(x_n, y_n,conn);

%radiation constant
h = 1;
% compute the element stiffness matrix and force vector
%and assemble them
[solution] = computeAndAssemble(elem, conn,ng,A,rhs,tnn,bc,x_n,y_n,h);

%sum length 
[anasol] = analyticalsolution(y_n,bc);

%errorNorm
%[err, nin] = errorNorm(solution,y);

%flag = 1;
%[err] = H1errornorm(ng,conn,elem,solution,flag,rhs,A,bc,max(x_n),max(y_n),sli)

toc

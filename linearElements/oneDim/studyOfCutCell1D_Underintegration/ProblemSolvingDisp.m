%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history
clc
clear all

%System length
L = 1;


%number of elements
ne = 64;

%number of nodes per element
nn = 2;

%create nodes
x_n1 = nodes(L,'true',ne,nn)

%weight
A1 = 3.0;

%BCs
%left BC
bc1 = 0.0;

%right BC
bc2 = 0.0;

%second variable 
A2 = 0.0;

%the total number of gauss points
ng = 2;

%constant RHS
K = 10;

%tolerance 
tol = 0.05;

% compute the element stiffness matrix and force vector
%and assemble them
[s,f] = computeAndAssemble(x_n1,A2,0,ng,K,L,nn);

%apply bcs
[s,f] = applyBCDisp(s,f,bc1, bc2);

%solve problem
tic
u_n1 = s\f;
toc

u_n1=[bc1;u_n1]


%System length
L = 1;


%number of elements
ne = 64;
elementLength = L/(ne)
cutElement = elementLength/2.0;
Lnew = L + cutElement
ne = ne+1;

%number of nodes per element
nn = 2;

%create nodes
x_n2 = nodes(Lnew,'true',ne,nn)

%weight
A1 = 3.0;

%BCs
%left BC
bc1 = 0.0;

%right BC
bc2 = 0.0;

%second variable 
A2 = 0.0;

%the total number of gauss points
ng = 2;

%constant RHS
K = 10;

%tolerance 
tol = 0.05;

% compute the element stiffness matrix and force vector
%and assemble them
[s,f] = computeAndAssembleCutCell(x_n2,A2,0,ng,K,L,nn);

%apply bcs
[s,f] = applyBCDispCutCell(s,f,bc1, bc2);

%solve problem
tic
u_n2 = s\f;
toc

u_n2=[bc1;u_n2]


plot(x_n1, u_n1, 'o',x_n2, u_n2,'+' ); %, x_n2, u_excat);

%plot(x_n, u_n, x_n2, u_n2);
%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history
clc
clear all

%System length
L = 1;


%number of elements
ne = 1000;
elementLength = L/(ne)
cutElement = elementLength/2.0;
Lnew = L + cutElement
ne = ne+1;

%number of nodes per element
nn = 2;

%create nodes
x_n = nodes(Lnew,'true',ne,nn)

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
[s,f] = computeAndAssembleCutCell(x_n,A2,0,ng,K,L,nn);

%apply bcs
[s,f] = applyBCDispCutCell(s,f,bc1, bc2);

%solve problem
tic
u_n = s\f;
toc

u_n=[bc1;u_n]

%evaluate the exact solution of this problem.
%x_n2 = 0:0.0001:1;

%u_excat = exactSolution(x_n2, K,L,A1,A2, bc1, bc2)

plot(x_n, u_n, '+'); %, x_n2, u_excat);

%plot(x_n, u_n, x_n2, u_n2);
    
    
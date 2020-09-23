%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history
clc
clear all
format long
%System length
L = 1.0; 
%number of volume elements
ne = 4;

%number of nodes per element
nn = 2;

%create vertices
x_n = vertices(L,'true',ne,nn)

%generate midpoints
c_n = centroids(x_n)

%weight
A = 3.0;

%BCs
%left BC
bc1 = 0.0;

%right BC
bc2 = 0.5;

%constant RHS
K = -3.0;

%length of cut cell
alpha = 0.9;

%length of volume element
DeltaX = x_n(2,1) - x_n(1,1);

%number of centroids
nc = length(c_n);

%get system matrix
s = SystemMatrix(c_n,alpha, DeltaX);
f = SystemRHS(DeltaX, alpha, K, A, nc, bc1, bc2);

%multiply with a minus sign, so the elements on the diagonal are postive.
%The way it is implemented the stiffness matrix is negative definite. We do care about definiteness whether postive or negative but it is nicer to have positive definiteness.

s = -1.0 * s;
f = -1.0 * f;

%solve problem
format long
u_n = s\f
format short

plot(c_n, u_n,'-o');


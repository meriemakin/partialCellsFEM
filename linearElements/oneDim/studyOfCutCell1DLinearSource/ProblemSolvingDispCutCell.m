%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history
clc
clear all
%format long

%System length
L = 1.0;

%number of elements
ne = 1024;

%h = [0.250000000000000,0.125,0.062500000000000,0.031250000000000,0.015625000000000,0.007812500000000,0.003906250000000,0.001953125000000,9.765625000000000e-04];
%dirichlet-dirichlet
%h0err = [0.015040107163804,0.003756601110705,9.394827965496672e-04,2.349085556766118e-04,5.873003838107202e-05,1.468270747913976e-05,3.670689757884591e-06,9.176732704281679e-07,2.294182967480236e-07];
%h1err = [0.041476939698894,0.020803922471080,0.010412634795572,0.005207807085728,0.002604099491843,0.001302074849272,6.510406007053279e-04,3.255206997472046e-04,1.627603999470796e-04];

%dirichlet-neumann
%h0err=[0.001541559124405,3.782741085182817e-04,9.411965760085815e-05,2.350181116744294e-05,5.873696117739108e-06,1.468314247785045e-06,3.670717227496756e-07,9.176749746685460e-08,2.294180634995586e-08];
%h1err= [0.041873446206524,0.020858463557356,0.010419765576442,0.005208718130046,0.002604214608408,0.001302089316267,6.510424139245294e-04,3.255209267025823e-04,1.627604283354380e-04];
%begin of system
xs = -0.5*(L/ne);
%xs = 0.0;

%end of system 
if(xs ~=0)
Le = xs + L + (L/ne);
else
Le = L + (L/ne);
end

%end of domain
alpha = 0.1;
% The way it is implemented here is that before the expansion to fitting purposes with the finite volume method, alpha is given with respect the last element. A new alpha has surely to be recalculated with respect to the new expanded domain. In order to maintian consistency of results, the old alpha has to be strictly between 0.5 and 1. But this means that alpha corresponding to the expanded domain can not exceed 0.5 which is kind of inconvenient mathematically, but physically it is consistent and correct. The reason why it can only be 0.5 at max. is because we expand the domain to half size of element length to the left and to the right. 

% If we use the expanded alpha as an actual alpha, i.e. that we ignore the fact that we have been expanding the domain, and that the last element totally belongs to the domain, then we can experiment with the effect of alpha onto the solution.

exp_alpha = alpha;

if(xs ~=0)
%0.2 is just specific for 0.9. I should think of a general formula to calculate the expanded alpha from  a given original alpha.
exp_alpha = 0.2;
end

Lc = L;
%Lc = (L - (L/ne)) + alpha * (L/ne);


%number of nodes per element
nn = 2;

if(xs~=0)
ne = ne + 1;
end

%Length of the element
elementL = (Le - xs)/ne;

%create nodes
x_n = nodes_vers2(xs,Le,ne,nn);
%plot(x_n);

%weight
A1 = 1.0;

%BCs

%left BC
bc1 = 0.0;

%right BC
bc2 = 0.5;

%second variable 
A2 = 0.0;

%the total number of gauss points
ng = 2;

%constant RHS
K = 0.0;

%tolerance 
tol = 0.05;

%compute the element stiffness matrix and force vector
%and assemble them
%first element to be used
if(xs~=0)
b = 1;
else
b = 1;
end

%1 if we want to scale the first and last element PDE contributions with the
%ratio of the physical domain to the entire element length. 
%0 if we want to consider the first and last element as they fully belong to the physical domain.
flag = 0;
[s,f] = computeAndAssembleCutCell(x_n,A1,0,ng,K,L,nn,b,exp_alpha,flag);
%add the constraint to the last element node
sold = s;

typeBC = 1;

if(typeBC == 1)
%#######################################################################
%                DIRICHLET BOUNDARY CONDITION ON CUT CELL
%#######################################################################

interD = 1;
if(interD == 1)
%Linear interpolation
[rows,cols] = size(s);
zerochunk = zeros(1,rows-2);
[weight1, weight2, rhs] = linearInterpolation(Le - elementL,Le,Lc,bc2);
%multiply with a minus sign, so the elements on the diagonal are postive.
%The way it is implemented the stiffness matrix is negative definite. We do care about definiteness whether postive or negative but it is nicer to have positive definiteness.
s(rows,:)=[zerochunk weight1 weight2];
f(rows)=rhs;

elseif(interD == 2)
%quadratic interpolation
[rows, cols] = size(s);
zerochunk = zeros(1,rows-3);
%position end of domain
[weight1, weight2, weight3, rhs] = quadraticInterpolation((Le -2*elementL),(Le - elementL),Lc,bc2,Le);
s(rows,:)=[zerochunk weight1 weight2 weight3];
f(rows)=rhs;

elseif(interD == 3)
%cubic interpolation
[rows, cols] = size(s);
zerochunk = zeros(1,rows-4);
elementL = L/ne;

 [weight1, weight2, weight3, weight4,rhs] = cubicInterpolation((L - 3*elementL),(L-2*elementL),(L-elementL),Lc,bc2,L);
s(rows,:)=[zerochunk weight1 weight2 weight3 weight4];
f(rows)=rhs;
end
%#######################################################################
%                NEUMANN BOUNDARY CONDITION ON CUT CELL
%#######################################################################

elseif (typeBC == 2)
interN = 2;
if(interN == 1)
%Derivative of quadratic interpolation
[rows,cols] = size(s);
zerochunk = zeros(1,rows-3);
[weight1, weight2, weight3, rhs] = quadInterpolationDerivative(Le - 2*elementL,Le - elementL,Lc,bc2);
s(rows,:)=[zerochunk weight1 weight2 weight3];
f(rows)=rhs;
elseif(interN == 2)
%Derivative of linear interpolation
[rows,cols] = size(s);
zerochunk = zeros(1,rows-2);
[weight1, weight2, rhs] = linearInterpolationDerivative(Le - elementL,Le,bc2);
s(rows,:)=[zerochunk weight1 weight2];
f(rows)=rhs;
end
end

if(xs ~=0)
%applyBCs
interBC = 1;
%linear interpolation
if(interBC == 1)
[rows,cols] = size(s);
zerochunk = zeros(1,rows-2);
[weight1, weight2, rhs] = linearInterpolation(xs,xs + elementL,0,bc1);
%multiply with a minus sign, so the elements on the diagonal are postive.
%The way it is implemented the stiffness matrix is negative definite. We do care about definiteness whether postive or negative but it is nicer to have positive definiteness.
weight1 = -1.0 * weight1;
weight2 = -1.0 * weight2;
rhs = -1.0 * rhs;
s(1,:)=[weight1 weight2 zerochunk];
f(1)=rhs;


%quadratic interpolation
elseif(interBC == 2)
[rows, cols] = size(s);
zerochunk = zeros(1,rows-3);
[weight1, weight2, weight3, rhs] = quadraticInterpolation(xs,xs + elementL,0,bc1,xs+2*elementL);
%multiply with a minus sign, so the elements on the diagonal are postive.
%The way it is implemented the stiffness matrix is negative definite. We do care about definiteness whether postive or negative but it is nicer to have positive definiteness.
weight1 = -1.0 * weight1;
weight2 = -1.0 * weight2;
weight3 = -1.0 * weight3;
rhs = -1.0 * rhs;
s(1,:)=[weight1 weight2 weight3 zerochunk];
f(1)=rhs;
end

else
[s,f] = applyBCDispCutCell(s,f,bc1, bc2);
end

%s
%f

%solve problem
format long
u_n = s\f;
format short


if(xs ~=0)
else
u_n = [bc1;u_n];
end

%figure
plot(x_n, u_n,'-o'); 

format long
max(svd(s))/min(svd(s));
%svd(s)
%svd(s-sold)


%########################################################################
%                         DERIVATIVE FEM SOLUTION
%########################################################################    

du_n = derivativeFEMSolution(nn,u_n, elementL);

%########################################################################
%                        ANALYTICAL SOLUTION
%########################################################################

u_nE = [];

if (typeBC == 1)
%analytical solution Dirichlet BC
u_nE = (-1/(6*A1))*x_n.^3+ ((bc2 - bc1 + (L^3/(6*A1)))/L)*x_n + bc1;  
elseif(typeBC == 2)
%analytical solution Neumann BC
u_nE = (-1/(6*A1))*x_n.^3 + ((bc2 + (L^2/(2)))/A1)*x_n + bc1;
end

format long
u_nE;  

%hold on
%plot(x_n, u_nE, '*'); 
%########################################################################
%                       DERIVATIVE ANALYTICAL SOLUTION
%########################################################################

%evaluated at element middle points (superconvergence points)
%analytical solution Dirichlet BC
du_nE = [];
if (typeBC == 1)
du_nE = -(1/(2*A1))*(x_n+(elementL/2.0)).^2 + ((bc2 - bc1 + (L^3/(6*A1)))/L);  
elseif (typeBC == 2)
%analytical solution Neumann BC
du_nE = -(1/(2*A1))*(x_n+(elementL/2.0)).^2 + ((bc2 + (L^2/(2)))/A1);
end

%figure
%plot(x_n, du_n,'o');
%hold on
%plot(x_n, du_nE,'*');

%########################################################################
%                       ERROR ANALYSIS
% 7.046872406736180e-05
% 7.046872406736180e-05
% 7.046872406736180e-05
% 7.046872406736180e-05
% 7.046872406736180e-05
% 7.046872406736180e-05
% 7.046872406736180e-05
%%########################################################################

%error made by the algorithm
[nin,err] = errorNorm(u_n, u_nE)
norm(u_n(2:end-1) - u_nE(2:end-1))/norm(u_nE(2:end-1));

%Dirichlet
%linear interpolation 
%nin  | err 
%-------------------------
% 4   | 0.004475171577105
% 8   | 0.001125432673009
% 16  | 2.817716264743138e-04
% 32  | 7.046872406736180e-05
% 64  | 1.761879423795508e-05
% 128 | 4.404799380047683e-06
% 256 | 1.101206146216950e-06
% 512 | 2.753019303877996e-07
%1024 | 6.882550722151122e-08


%quadratic interpolation
%nin | err
% 4  | 5.459150335694010e-04
% 8  | 6.992455589153291e-05
% 16 | 8.792442759400491e-06
% 32 | 1.100670403711856e-06
% 64 | 1.376342514322420e-07
% 128| 1.720565624686261e-08
% 256| 2.151171677959698e-09
% 512| 2.657501369947409e-10
%1024| 4.136579276853673e-11

%Neumann
%derivative of quadratic interpolation
%nin  | err
%4    | 7.973599423124790e-04
%8    | 2.759737135061665e-04
%16   | 8.106169357410469e-05
%32   | 2.185835957662591e-05
%64   | 5.667702143601524e-06 
%128  | 1.442533624155729e-06 
%256  | 3.638471396743914e-07 
%512  | 9.136379799487464e-08 
%1024 | 2.289115150959034e-08 

%derivative of linear interpolation
%nin  | err
%4    | 0.001491723859035 
%8    | 3.751442243362687e-04 
%16   | 9.392387549209402e-05 
%32   | 2.348957469170692e-05 
%64   | 5.872931411391047e-06 
%128  | 1.468266462601664e-06 
%256  | 3.670686948194450e-07 
%512  | 9.176730944623162e-08 
%1024 | 2.294187678539610e-08


flagerror = 1;
[errH0H1] = errornorms(ng,x_n,nn,u_n,flagerror,K,A1,bc1,bc2,L,typeBC) 
h = x_n(2) - x_n(1)
%Results
%h=[];
%dirichlet-dirichlet-problem
%H0-error

%H1-error

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%dirichlet-neumann problem
%H0-error

%H1-error

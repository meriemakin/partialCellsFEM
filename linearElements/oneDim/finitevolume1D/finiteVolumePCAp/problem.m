%clear command window and variable history
clc
clear all
format long
%System length
L = 1.0;

%number of volume elements
ne = 512;

%number of nodes per element
nn = 2;

%create vertices
x_n = vertices(L,'true',ne,nn);

%generate midpoints
c_n = centroids(x_n);

%weight
A = 3.0;

%BCs
%left BC
bc1 = 0.0;

%right BC
bc2 = 0.5;
% 0 : dirichlet, 1: neumann
flag = 1;

%constant RHS
K = 1;

%length of cut cell
alpha = 0.9;

%length of volume element
DeltaX = x_n(2,1) - x_n(1,1);

%number of centroids
nc = length(c_n);

%get system matrix
s = SystemMatrix(c_n,alpha, DeltaX,flag);
f = SystemRHS(DeltaX, alpha, K, A, nc, bc1, bc2,flag);


%multiply with a minus sign, so the elements on the diagonal are postive.
%The way it is implemented the stiffness matrix is negative definite. We do care about definiteness whether postive or negative but it is nicer to have positive definiteness.

s = -1.0 * s;
f = -1.0 * f;

%solve problem
format long
u_n = s\f;

plot(c_n, u_n,'o');


%analytical solution
if(flag == 0)
	c_1 = (1/L)*(bc2 - bc1 + (K/(2*A))*L^2);
	c_2 = bc1;
elseif(flag == 1)
	c_1 = bc2 + ((K/A)*L);
	c_2 = bc1;
end
u_ex = [];
u_exder = [];
for i=1:length(c_n)
      	u_ex = [u_ex; -(K/(2*A)) * (c_n(i,1))^2 + c_1*c_n(i,1) + c_2];
      	u_exder = [u_exder; -(K/(A)) * c_n(i,1) + c_1];
end 
hold on
plot(c_n,u_ex,'x')
h = L/ne;

%calculate h0-error
h0err = 0.0;
for i=1:length(u_ex)
	%midpoint rule
	h0err = h0err + ((u_n(i,1) - u_ex(i,1))*h)^2.0;
end
h0err = sqrt(h0err);
h0err
%results
%dirichlet-dirichlet
%h = [0.250000000000000, 0.125000000000000, 0.062500000000000, 0.031250000000000, 0.015625000000000, 0.007812500000000, 0.003906250000000, 0.001953125000000]
%h0err = [0.003671935653009, 0.001289399209890, 4.536700401095130e-04, 1.599523716399855e-04, 5.646837339486031e-05, 1.994943196637091e-05, 7.050473425775824e-06, 2.492235394599598e-06];
%dirichlet-neumann
%h0err=[0.002386758174456,8.488544798439725e-04,3.005564015727554e-04,1.063016803192136e-04,3.758676104094064e-05,1.328923098403537e-05,4.698479564306181e-06,1.661165773097870e-06];
%poly=polyfit(log(h),log(h0err),1);
%slope = poly(1);

%calculate h1-error
h1err = 0.0;
for i=1:length(u_ex)
        %midpoint rule
        h1err = h1err + ((u_n(i,1) - u_ex(i,1))*h)^2.0 + (u_exder(i,1)*h)^2.0;
end
h1err = sqrt(h1err);
h
h1err

%results
%dirichlet-neumann
%h = [0.250000000000000, 0.125000000000000, 0.062500000000000, 0.031250000000000, 0.015625000000000, 0.007812500000000, 0.003906250000000, 0.001953125000000];
%h1err = [0.250026964768682,0.176781397636523,0.125000823263310,0.088388492377465,0.062500025509412,0.044194178326787,0.031250000795347,0.022097087052624];
%dirichlet-neumann
%h1err=[0.336581263017690,0.238108376640030,0.168387384486753,0.119071278949292,0.084196712892274,0.059536173429521,0.042098450835790,0.029768103400602];

%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history
clc
clear all

%System length
L = 1;

%number of elements
ne = 4;

%begin of system
xs = -0.5*(L/ne);
%xs = 0.0;

%end of domain
%alpha is calculated with respect of the non-expanded domain. In order to match the finite volume method, the domain is expanded to the right and to the left Hence, alpha changes. It has to be thought of a general formula.
alpha = 0.9;
exp_alpha = alpha;
if(xs ~=0)
%0.2 is just specific for 0.9. I should think of a general formula to calculate the expanded alpha from  a given original alpha.
exp_alpha = 0.2;
end
%alpha factor changes but the end of the domain should remain consistent.
Lc = (L - (L/ne)) + alpha * (L/ne);

%end of system
if(xs ~=0) 
Le = xs + L + (L/ne);
else
Le = L;
end 

%number of nodes per element
nn = 2;

if(xs~=0)
ne = ne + 1;     
end

%Length of the element
elementL = (Le - xs)/ne;


%create nodes
%x_n = nodes(L,'true',ne,nn)
x_n = nodes_vers2(xs,Le,ne,nn)

%weight
A1 = 3.0;

%BCs
%left BC
bc1 = 0.0;

%right BC
bc2 = 0.5;

%second variable 
A2 = 0.0;

%the total number of gauss points
ng = 2;

%constant RHS %body force
K = 0;

%tolerance 
tol = 0.05;

%penalty coefficient 
p = 1000;


%kappa = 0.5
%u_n =
%
%         0
%    0.0620
%    0.1240
%    0.1860
%    0.2480
%    0.3100
%    0.3720
%    0.4340
%    0.4960
%    0.5580
%    0.6199
%    0.6819
%    0.7435

%kappa = 1.0
%u_n =
%
%         0
%    0.0640
%    0.1279
%    0.1919
%    0.2558
%    0.3198
%    0.3837
%    0.4477
%    0.5116
%    0.5756
%    0.6395
%    0.7035
%    0.7353

% compute the element stiffness matrix and force vector
%and assemble them
[s,f] = computeAndAssembleCutCell(x_n,A2,0,ng,K,Lc,nn,bc2,p,exp_alpha);

%multiply with a minus sign, so the elements on the diagonal are postive.
%The way it is implemented the stiffness matrix is negative definite. We do care about definiteness whether postive or negative but it is nicer to have positive definiteness.

s = -1.0 * s;
f = -1.0 * f;

%apply bcs
if (xs ~=0)
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
else
[s,f] = applyBCDispCutCell(s,f,bc1, bc2);
end

%solve problem
format long
u_n = s\f
format short

if(xs==0)
u_n=[bc1;u_n]
end
%evaluate the exact solution of this problem.
%x_n2 = 0:0.0001:1;

%u_excat = exactSolution(x_n2, K,L,A1,A2, bc1, bc2)

plot(x_n, u_n, '+-'); %, x_n2, u_excat);

%plot(x_n, u_n, x_n2, u_n2);
    
pen=[10000,1000,100,10,5,1];
tipDeflection = [0.5003,0.5016,0.5156,0.7145,1.2508,-0.2500];
%plot(log(pen),log(tipDeflection),'+-');

numberElem =[4,10,100,500,1000];
tipDefltectionPen1000=[0.5255,0.5127, 0.5027,0.5018,0.5016];  
tipDefltectionPen100=[0.5404,0.5272,0.5168,0.5157,0.5156]; 
tipDefltectionPen10=[0.7545,0.7353,0.7168,0.7148,0.7145]; 
tipDefltectionPen1=[-0.2547,-0.2495,-0.2497,-0.2499,-0.2500]; 
%plot(numberElem,tipDefltectionPen1000,'+-',numberElem,tipDefltectionPen100,'o-',numberElem,tipDefltectionPen10,'x-',numberElem,tipDefltectionPen1,'.-');

%Element Stiffness
%temp =

%    0.2500    0.2500
%    0.2500    0.2500

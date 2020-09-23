%clear command window and variable history
clc
clear all
%format long

%System length
L = 1.0;

%number of elements
ne = 128; 


%begin of system
xs = -0.5*(L/ne);
%xs = 0.0;

% The way it is implemented here is that before the expansion to fitting purposes with the finite volume method, alpha is given with respect the last element. A new alpha has surely to be recalculated with respect to the new expanded domain. In order to maintian consistency of results, the old alpha has to be strictly between 0.5 and 1. But this means that alpha corresponding to the expanded domain can not exceed 0.5 which is kind of inconvenient mathematically, but physically it is consistent and correct. The reason why it can only be 0.5 at max. is because we expand the domain to half size of element length to the left and to the right. 

% If we use the expanded alpha as an actual alpha, i.e. that we ignore the fact that we have been expanding the domain, and that the last element totally belongs to the domain, then we can experiment with the effect of alpha onto the solution.


factor = 0.5;
Lc = L;
Le = Lc + factor *(L/ne);

%number of nodes per element
nn = 2;

if(xs~=0)
ne = ne + 1;
end

%Length of the element
elementL = (Le - xs)/ne;

%create nodes
x_n = nodes_vers2(xs,Le,ne,nn);
exp_alpha = (x_n(end)-Lc)/(x_n(end)-x_n(end-1));

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
K = 5.0;

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

typeBC = 2;

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
%plot(x_n, u_n,'-o'); 

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
u_nE = (-K/A1)*0.5 * (x_n.^2) + (((bc2 - bc1)/Lc) + (Lc/2)*(K/A1))*x_n + bc1;  
elseif(typeBC == 2)
%analytical solution Neumann BC
u_nE = 0.5 * (-K/A1) * x_n.^2 + (bc2 + (K/A1)* Lc) * x_n + bc1;
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
du_nE = (-K/A1)*(x_n+(elementL/2.0)) + (((bc2 - bc1)/Lc) + (Lc/2)*(K/A1));  
elseif (typeBC == 2)
%analytical solution Neumann BC
du_nE =  (-K/A1)*(x_n+(elementL/2.0)) + (bc2 + (K/A1)*Lc);
end

%figure
%plot(x_n, du_n,'o');
%hold on
%plot(x_n, du_nE,'*');

%########################################################################
%                       ERROR ANALYSIS
%########################################################################

%error made by the algorithm
[nin,err] = errorNorm(u_n, u_nE);

norm(u_n(2:end-1) - u_nE(2:end-1))/norm(u_nE(2:end-1));

%error results for linear interpolation on both sides
% err = [0.031250000000000,0.007812500000000,0.001953125000000,4.882812500000000e-04,1.220703125001110e-04,3.051757812500000e-05,7.629394531250000e-06, 1.907348632701478e-06,4.768371582031250e-07,1.192092895507812e-07,2.980232238769531e-08,7.450580596923828e-09];
% nin = [2,4,8,16,32,64,128,256,512,1024,2048,4096];


%error results for quadratic interpolation on both sides
%err = [4.163336342344337e-17,1.296705798292663e-16,1.112825109839122e-15,2.465150479580291e-15,2.081267016368349e-14,4.223409274216328e-14,3.404037170412835e-13,6.676058021833596e-13,5.472485393445345e-12,1.088767981445528e-11, 8.732092397230042e-11, 1.748866421536470e-10, 1.396787000525363e-09];
%nin = [2,4,8,16,32,64,128,256,512,1024,2048,4096,8192] ;
%conditionnumberafter = [4.887839941639816, 9.611477771805005, 31.047976669642463,1.137030096500675e+02,4.348998553018432e+02,1.699906386269271e+03,6.720031545111005e+03,2.672058006095946e+04,1.065627977885945e+05,4.256116804910579e+05,1.701167217338147e+06,6.802109365444869e+06,2.720331793118353e+07];

%########################################################################
%########################################################################
%              Eigenvalue Analysis- condition number
%########################################################################
%########################################################################
%h = [0.5,0.25, 0.1250,0.0625,0.03125,0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625000000000e-04];
%before applying boundary conditions ignoring the singular values that are close to the machine epsilon (i dont know if this makes sense)
%conditionnumber = [5.828427124746193,13.928203230275521,39.863458189061546,1.306460956438605e+02,4.678426288390451e+02];
%taking into account all singular values
%conditionnumber = [4.085925719606174e+16,1.288222685639135e+16,1.102518255522341e+16,1.610545333263965e+17,4.216490311080406e+16, 2.842184018423543e+16, 2.849914728689063e+16,2.069591063438704e+16,9.754538101110034e+16,1.384525765787646e+16];
%after applying boundary conditions
%conditionnumberafter = [4.887839941639816, 9.611477771805005, 31.047976669642463,1.137030096500675e+02,4.348998553018432e+02,1.699906386269271e+03,6.720031545111005e+03,2.672058006095946e+04,1.065627977885945e+05,4.256116804910579e+05,1.701167217338147e+06,6.802109365444869e+06,2.720331793118353e+07];

%########################################################################
%                    MORE ACCURATE ERROR ANALYSIS
%########################################################################

%flagerror
%if 0 H0 norm
%if 1 H1 norm
flagerror = 0;
[errH0H1] = errornorms(ng,x_n,nn,u_n,flagerror,K,A1,bc1,bc2,Lc,typeBC)
h = x_n(2)- x_n(1)
%dirichlet-dirichlet problem with linear interpolation on cut cell side and 
%quadratic interpolation on regular end

%alpha 0.3
%h = [0.260000000000000,0.127777777777778,0.063235294117647, 0.031346153846154, 0.015673076923077, 0.007824612403101, 0.003909289883268, 0.001953886452242, 9.767530487804878e-04];
%h0err=[0.012813424890144,0.002880591387497,6.759676819054190e-04,1.631861495490018e-04,4.005147342894875e-05,9.918469891827951e-06,2.467736734695221e-06,6.154420052126405e-07,1.536745480288784e-07];
%h1err=[0.406763192576320,0.192320315751822,0.093245655459947,0.045872239823587,0.022745481357153,0.011324690217479,0.005650282534716,0.002822119308197,0.001410303370373];

%alpha 0.5
%h=[0.235000000000000,0.120833333333333,0.061397058823529,0.030965909090909, 0.015552884615385,0.007794331395349, 0.003901690175097,0.001951982821637,9.762766768292683e-04];
%h0err=[0.013810679320050,0.003355391947931,8.264199223465982e-04,2.050338172479952e-04,5.106092581495934e-05,1.274046861480681e-05,3.182017282260548e-06,7.951165533834953e-07,1.987306497155404e-07];
%h1err= [0.382981866232960, 0.186004965516491, 0.091613410340752, 0.045456972964329, 0.022640727733147, 0.011298382327285, 0.005643690474954, 0.002820469397763,0.001410234698881];


%alpha 0.8 
%h=[0.235000000000000,0.120833333333333,0.061397058823529,0.030965909090909, 0.015552884615385,0.007794331395349, 0.003901690175097,0.001951982821637,9.762766768292683e-04];
%h0err=[0.007855905720679,0.001892775287726,4.631917518645474e-04,1.144608492829513e-04,2.844193367674319e-05,7.088414361767122e-06,1.769318498946266e-06,4.419800260025640e-07,1.104597563844063e-07];
%h1err=[0.350371053978160,0.181917291828379,0.090541939896552,0.045182310025024,0.022571170948569,0.011280879006387,0.005639300202984,0.002819370013392,0.001409615581383];



%dirichlet-neumann problem with linear interpolation on cut cell side and 
%quadratic interpolation on regular end
%alpha 0.3
%h0err = [0.150894951008243,0.073522775364362,0.036382429600508,0.018111295080325, 0.009037672459825, 0.004514603501619, 0.002256276884036, 0.001127886431407, 5.638807438995596e-04];
%h = [0.260000000000000,0.127777777777778,0.063235294117647, 0.031346153846154, 0.015673076923077, 0.007824612403101, 0.003909289883268, 0.001953886452242, 9.767530487804878e-04];
%h1err=[0.474153764356928,0.233890378601813,0.116211387348007,0.057931219576635,0.028923191738706,0.014451146566361,0.007222980690159,0.003610844685658,0.001805261242321];

%alpha 0.5

%h=[0.250000000000000, 0.125000000000000, 0.062500000000000, 0.031250000000000, 0.015625000000000, 0.007812500000000, 0.003906250000000, 0.001953125000000,9.765625000000000e-04];
%h0err=[0.013810679320048, 0.003355391947931, 8.264199223492949e-04, 2.050338172587852e-04, 5.106092580960570e-05,1.274046862548549e-05,3.182017196902461e-06,7.951165507143328e-07,1.987308202893025e-07];
%h1err=[0.382981866232960,0.186004965516491,0.091613410340752,0.045456972964329,0.022640727733147,0.011298382327285,0.005643690474954,0.002820469397763,9.765625000000000e-04] ;



%alpha 0.8
%h=[0.235000000000000,0.120833333333333,0.061397058823529,0.030965909090909, 0.015552884615385,0.007794331395349, 0.003901690175097,0.001951982821637,9.762766768292683e-04];
%h0err=[0.165843781083111,0.103476060734010,0.052842145033830,0.026730299950875,0.013446870450097,0.006744435533654,0.003377540341641,0.001690109937583,8.453909269971209e-04];
%h1err=[0.544545557776108,0.279086670763414,0.140111039687780,0.070240277065380,0.035172237287318,0.017599920762308,0.008803510651104,0.004402655508781,0.002201554203943];




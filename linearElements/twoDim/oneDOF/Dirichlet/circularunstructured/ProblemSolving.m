%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history

tic

clc
clear all
format long

%Radius of the Circle
R = 2;

%center coordinates
xC = 1;
yC = 1;
%refinement level 
rl = 7; 

%Errors after neighbor averaging
%h  = [1, 0.5,0.2500,0.1250,0.062500000000000,0.031250000000000,0.01562500000000, 0.007812500000000];
%H0 = [0.475128864304853, 0.135003169636846, 0.035019892203699, 0.008752253863575, 0.002158901176438, 5.366181723234731e-04];
%H1 = [];
%Energy = [];

%number of gauss points in one direction
ng = 2;


%Material factor in the DE
A = 1.0;

%The rhs of the system
rhs = 1.0;

%BCs on the boundary of the circle 
bc = 0.0;

%create elements
[elem,conn,h,flags,nn,x_n,y_n,nel] = elements(xC,yC,R,rl);
%figure
%for i=1:nel
    %elem(1,:,i)
    %elem(2,:,i)
%    for j=1:4
%        plot(elem(1,j,i),elem(2,j,i),'+');
%        hold on
%    end
%end

ncycles = 40;
[x_n,y_n,elem] = neighborAveraging(flags, elem,conn,x_n,y_n,nn,ncycles);

%figure
%plot(x_n,y_n,'+');
%figure
%for i=1:nel
    %elem(1,:,i)
    %elem(2,:,i)
%    for j=1:4
%        plot(elem(1,j,i),elem(2,j,i),'+');
%        hold on
%    end
%end

h
% compute the element stiffness matrix and force vector
%and assemble them
[solution] = computeAndAssemble(elem, conn,ng,A,rhs,nn,bc,flags);

%calculate the analytical solution
%[y] = analyticalsolution(x_n,y_n,bc);
%series last index

[anasol] = analyticalsolution(rhs,A,x_n,y_n,xC,yC,R,bc); 

%errorNorm
%[err, nin] = errorNorm(solution,y);

flag = 1;
[err] = H1errornorm(ng,conn,elem,solution,flag,rhs,A,bc,xC,yC,R)

%Errors

%h  = [1, 0.5,0.2500,0.1250,0.062500000000000,0.031250000000000,0.01562500000000, 0.007812500000000];
%H0 = [0.479188979127252 ,0.1333, 0.0344,0.008667743112118,0.002171773278330,5.432574322855775e-04,1.358345247342379e-04,3.395991913712670e-05];
%H1 = [0.857449795402268, 0.415817544458165, 0.205922549005989,0.102747961793460,0.051353345985281,0.025674711991492,0.012837169951385];
%Energy = [0.711054199001998, 0.393858034431364, 0.203032317991283, 0.102381706774466, 0.051307402435823, 0.025668963890432, 0.012836451275280];


%Errors after neighbor averaging
%H0_2 = [0.475128864304853, 0.135003169636846, 0.035019892203698,0.008752253863575,0.002158901176432,5.366181723205258e-04 ];
%H1_2 = [0.847661195532879,0.420396992992975,0.208994939352206,0.103203157923861, 0.050710680248005,0.025203671187034 ];
%Energy_2 = [0.701984376404923, 0.398130350394868,0.206040024812835,0.102831366118344, 0.050664704054458,0.025197957898236 ];

toc

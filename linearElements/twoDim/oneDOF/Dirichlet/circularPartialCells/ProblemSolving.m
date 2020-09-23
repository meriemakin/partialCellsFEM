%---------------------------------------------------------------------------%
%                           This is a script file.                          %             
%---------------------------------------------------------------------------%
%clear command window and variable history

tic

clc
clear all

%Length
Lx = 2;

%width
Ly = Lx;

%Number of elements in x-direction
Nx = 512; 
%Number of elements in y-direction
Ny = Nx;


%number of gauss points in one direction
ng = 2;

%coordinate of center of circular physical domain
xC = Lx / 2.0;
yC = Ly / 2.0;

%radius of circular physical domain
radius = (Lx/2.0);

%Material factor in the DE
A = 1.0;

%The rhs of the system
rhs = 1.0;

%BCs on the boundary of the circular physical domain 
bc = 0.0;

%create nodes
[x_n,y_n] = nodes(Lx, Ly, Nx, Ny);
%plot(x_n, y_n, 'o');
%hold on;
%Center = [xC,yC];
%circle(Center,radius,1000,'b-');
%drawnow;
%pause

%total number of nodes
tnn = length(x_n);

%the connectivity
conn = connectivity(Nx, Ny,2,2);

%create elements
elem = elements(x_n, y_n,conn);

%flag elements whether they are inside of physical domain, outside of physical domain or partially covered by physical domain.
flagsE = elementflags(elem, xC, yC, radius);

%bias the way the collocation method works
%0 dont bias
%1 bias
bias  = 0;

% compute the element stiffness matrix and force vector
%and assemble them
[solution] = computeAndAssemble(elem, conn,ng,A,rhs,tnn,flagsE,xC,yC,radius,bc,x_n,y_n,bias);

%calculate the analytical solution
[y] = analyticalsolution(rhs,A,x_n,y_n,xC,yC,radius,bc);

%errorNorm
%[err, nin] = errorNorm(solution,y);

flag = 1; 
err =  H1errornorm(ng,conn,elem,solution,flag,rhs,A,bc,xC,yC,radius,flagsE)
h = x_n(2)-x_n(1)

%some results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolation techniques involves only one direction since the intersection
%point lies between two nodes. So a linear onedimensional interpolation is 
%done. For this algorithm, here is the error analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%using the function errornorm
%nin = [5, 13, 49, 197,797,848, 1257, 3209,12853,51433,86584];
%err = [0.022613849221600, 0.012440258770754, 0.003636812889227, 7.512968349467671e-04,2.211819084109193e-04, 1.527300058654243e-04, 1.160100600789700e-04, 4.266346728977254e-05,1.173832783656672e-05, 3.024078476989026e-06, 1.847789491841085e-07];
%sqrterr = [0.140054426101652,0.108713340935834,0.059679142435583, 0.027413343970278,0.014915513311829,0.006615323826820,0.003463160780035,0.001750877985589 ];

%figure
%loglog(nin, err, '+-');

%using the function H1errornorm
%H0 Norm
%h = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625];
%err = [0.050363188827124,0.016305493915506,0.003767726874018,0.001077008500860,2.054219166814060e-04,6.832911361378792e-05, 1.522770561721046e-05];


%H1 Norm
%err = [0.302058508577937,0.186362949908247,0.093417446498534,0.046104752983967,0.022811857319487,0.011353995233651,0.005661145545804,0.002828739200443];
%h = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125];


%Energy
%err= [0.251695319750813,   0.170057455992741 ,  0.089649719624516, 0.045027744483107,   0.022606435402806,   0.011285666120037, 0.00564591784018];
%h = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Error analysis for a bilinear interpolation scheme. All nodes are involved
%the interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%   NO BIASING
%%%%%%%%%%%%%%%%%

%H0error = [0.050363188827124, 0.012832476965068, 0.001573908860193, 3.978374469148100e-04, 5.387958770386010e-05, 1.615582208101466e-05,2.584625372372554e-06, 5.514063362883122e-07,1.440695550443632e-07];
%H1error = [0.302058508577937, 0.186951492965866, 0.094031864634782, 0.045606951190493,0.022670808488914,0.011318605361119,0.005658257049788,0.002827370282523,0.001413706837007 ];
%En = [0.251658508577937,   0.174151492965866,   0.092431864634782, 0.045209111190493,   0.022616928488914,   0.011302449361119,  0.005655672449788,   0.002826818876187];
%h = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125,0.003906250000000];


%%%%%%%%%%%%%%%%%%
%    WITH BIASING
%%%%%%%%%%%%%%%%%%
%h = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125,0.003906250000000];
%H0error_1 = [0.050363188827124, 0.011608049386764,0.001641280872339,3.775202179678777e-04,4.797418632858967e-05,1.178159085166205e-05,1.694223450582157e-06,2.974091297341692e-07, 5.716774911249880e-08];
%H1error_1 = [0.302058508577937,0.186587852385474,0.092540898513000,0.045293566270756,0.022600411568611,0.011299175517420,0.005652306010347,0.002826240036083,0.001413273261189];


toc

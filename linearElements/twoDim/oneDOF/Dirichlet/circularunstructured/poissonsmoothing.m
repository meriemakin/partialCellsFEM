%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%                      Poisson Smoothing
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x_n,y_n, elem] = poissonsmoothing(flags, elem, conn, x_n, y_n, tnn, ncycles);

%create nodes
x = [0, 4, 8, 12, 16;
     0, 3, 9, 11, 16;
     0, 5, 7, 13, 16;
     0, 4, 8, 12, 16];

y = [12, 12, 12, 12, 12;
     8, 9, 7, 8, 8;
     4, 5, 4, 3, 4;
     0, 0, 0, 0, 0];


%optimal location:
x_2 = [0, 4, 8, 12, 16;
       0, 4, 8, 12, 16;
       0, 4, 8, 12, 16;
       0, 4, 8, 12, 16];

y_2 = [12, 12, 12, 12, 12;
       8, 8, 8, 8, 8;
       4, 4, 4, 4, 4;
       0, 0, 0, 0, 0];

indices = [4, 8, 12 , 16, 20;
           3,7,11,15,19;
           2,6,10,14,18;
           1,5,9,13,17]
  
%calculate the distance of the two nodes from their placement as a function of the number of cycles
d1 = sqrt((x(3,2)-x_2(3,2))^2+(y(3,2)-y_2(3,2))^2);
d2 = sqrt((x(2,2)-x_2(2,2))^2+(y(2,2)-y_2(2,2))^2);
d3 = sqrt((x(3,3)-x_2(3,3))^2+(y(3,3)-y_2(3,3))^2);
d4 = sqrt((x(2,3)-x_2(2,3))^2+(y(2,3)-y_2(2,3))^2);
d5 = sqrt((x(3,4)-x_2(3,4))^2+(y(3,4)-y_2(3,4))^2);
d6 = sqrt((x(2,4)-x_2(2,4))^2+(y(2,4)-y_2(2,4))^2);

distance = [0,0,0,0,0,d1,d2,0,0,d3,d4,0,0,d5,d6,0,0,0,0,0];

%poisson averaging
num_cycles = 20;

for k=1:num_cycles
        dist_temp = zeros(1,20);
    for i=2:3
        for j=2:4
        x_xi = 0.5*(x(i+1,j)-x(i-1,j));
        x_eta = 0.5*(x(i,j+1)-x(i,j-1));
        y_xi = 0.5*(y(i+1,j)-y(i-1,j));
        y_eta = 0.5*(y(i,j+1)-y(i,j-1)); 
        alpha = x_eta^2 + y_eta^2;
        beta = x_xi*x_eta+y_xi*y_eta;
        gamma = x_xi^2+y_xi^2;
        x_xieta = 0.25*(x(i+1,j+1)+x(i-1,j-1)-x(i+1,j-1)-x(i-1,j+1)); 
        y_xieta = 0.25*(y(i+1,j+1)+y(i-1,j-1)-y(i+1,j-1)-y(i-1,j+1));
        x(i,j) = -2*beta*x_xieta-alpha*(x(i+1,j)+x(i-1,j))-gamma*(x(i,j+1)+x(i,j-1));
        x(i,j) = x(i,j)/(-2*alpha - 2*gamma);
        y(i,j) = -2*beta*y_xieta-alpha*(y(i+1,j)+y(i-1,j))-gamma*(y(i,j+1)+y(i,j-1));
        y(i,j) = y(i,j)/(-2*alpha - 2*gamma);
        %calculate the distance of the two nodes from their placement as a function of the number of cycles
        dist_temp(indices(i,j)) = sqrt((x(i,j)-x_2(i,j))^2+(y(i,j)-y_2(i,j))^2);
        end
    end
    distance = [distance; dist_temp];
end
% plot the distances for every smoothing cylce
node_indices = 1:20;

for i=1:num_cycles
    plot(node_indices, distance(i,:), 'o-');
    hold on
end

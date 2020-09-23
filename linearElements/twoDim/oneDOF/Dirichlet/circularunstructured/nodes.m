%---------------------------------------------------------------------------%                          
%                           generating nodes                                %
%---------------------------------------------------------------------------%

%R is the radius of the circle
%Nx is the number of elements in x-direction
%Ny is the number of elements in y-direction
function [x_n,y_n] = nodes(R,xC,yC, Nr, Ntheta);

x_n=zeros((Nr+1)*(Ntheta+1),1);
y_n=zeros((Nr+1)*(Ntheta+1),1);
theta = (2*pi)/Ntheta;
step = R/Nr;

for j=0:Ntheta
    for i= 1:(Nr+1)
        x_n((i+ j*(Nr+1)),1) = (xC + step * (i-1))*cos(j*theta);
        y_n((i+ j*(Nr+1)),1) = (yC + step * (i-1))*sin(j*theta);
    end
end

plot(x_n, y_n, '+');



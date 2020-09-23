%---------------------------------------------------------------------------%                          
%                           generating nodes                                % 
%---------------------------------------------------------------------------%

%Lx inner radius of the arch
%Ly is the outer radius of the arch
%Nx is the number of elements in x-direction
%Ny is the number of elements in y-direction
function [x_n,y_n] = nodes(Lx, Ly, Nx, Ny);

x_n=zeros((Nx+1)*(Ny+1),1);
y_n=zeros((Nx+1)*(Ny+1),1);

stepX = Lx/Nx;
stepY = Ly/Ny;

for j=1:(Ny+1)
    	for i= 1:(Nx+1)    
             y_n((i+ (j-1)*(Nx+1)),1) = stepY * (j-1);	
             x_n((i+ (j-1)*(Nx+1)),1) = stepX * (i-1);
    	end
end

%plot(x_n, y_n, '+');



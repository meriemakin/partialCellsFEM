%---------------------------------------------------------------------------%                          
%                           Error norm                                      %
%---------------------------------------------------------------------------%

function [err] = H1errornorm(ng,conn,elem,solution,flag,f,A,bc,xC,yC,radius);


%flag
%if 0 H0 norm
%if 1 H1 norm
%if 2 energy norm

[nelements, numberofnodes] = size(conn);

err  = 0;

for k =1:nelements
    
    	el = elem(:,:,k);
    	connEl = conn(k,:);
    	u_nE = elementSolution(solution,connEl);

    	for i =1:ng
        	for j=1:ng
            		[xi_g,eta_g,w1,w2] = gaussianquadrature(ng,i,j);
            		phi = shapefunction(xi_g, eta_g);
            		[x,y] = mapping(phi, el);
            
                        Dphi = shapefunctionFirstDerivative(xi_g,eta_g);
                        F = deformationGradient(el,Dphi);
                	jacob = jacobian(F);
                
                	if(flag == 0 || flag == 1)
                           u = analyticalsolutionSinglePoint(f,A,bc,x,y,xC,yC,radius);
                   	   err = err + (phi*u_nE - u)^2*jacob*w1*w2; 
                	end
 
                	if(flag == 1)
                   	   mat = u_nE'*Dphi*F^-1;
                   	   [dudx,dudy] = derivativeAnalyticalsolution(f,A,x,y,xC,yC);
                   	   err = err + ((mat(:,1) - dudx)^2 + (mat(:,2) - dudy)^2)*jacob*w1*w2; 
                	end

                        if(flag == 2)
                           mat = u_nE'*Dphi*F^-1;
                           [dudx,dudy] = derivativeAnalyticalsolution(f,A,x,y,xC,yC);
                           err = err + A*((mat(:,1) - dudx)^2 + (mat(:,2) - dudy)^2)*jacob*w1*w2;
                        end
         	end
    	end
    end

%take square root of everything
err = sqrt(abs(err));


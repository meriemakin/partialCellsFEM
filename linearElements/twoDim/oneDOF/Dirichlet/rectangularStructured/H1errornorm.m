%---------------------------------------------------------------------------%                          
%                           Error norm                                      %
%---------------------------------------------------------------------------%

function [err] = H1errornorm(ng,conn,elem,solution,flag,f,A,bc,a,b,n);


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
        %u_nE = [u_nE(1); u_nE(3); u_nE(4) ;u_nE(2)];
    	for i =1:ng
        	for j=1:ng
            		[xi_g,eta_g,w1,w2] = gaussianquadrature(ng,i,j);
            		phi = shapefunction(xi_g, eta_g);
            		[x,y] = mapping(phi, el);
                		
                        Dphi = shapefunctionFirstDerivative(xi_g,eta_g);
                	F = deformationGradient(el,Dphi);
                	jacob = jacobian(F);
                        %jacob = abs(jacob);        
                	
                        if(flag == 0 || flag == 1)
                           u = analyticalsolutionSinglePoint(x,y,a,b,f,A,n);
                   	   err = err + (phi*u_nE - u)^2*jacob*w1*w2; 
                	end
 
                	if(flag == 1)
                   		mat = (Dphi*F^-1)'*u_nE;
                                %mat = u_nE'*Dphi*F^-1;
                   		[dudx,dudy] = derivativeAnalyticalsolution(x,y,a,b,f,A,n);
                                err = err + (((mat(1) - dudx)^2 + (mat(2) - dudy)^2)*jacob*w1*w2); 
                	end

                        if(flag == 2)
                               mat = u_nE'*Dphi*F^-1;
                               [dudx,dudy] = derivativeAnalyticalsolution(x,y,a,b,f,A,n);
                               err = err + A*((mat(:,1) - dudx)^2 + (mat(:,2) - dudy)^2)*jacob*w1*w2;
                        end
      		end 
        
       	end
end

%take square root of everything
err = sqrt(abs(err));


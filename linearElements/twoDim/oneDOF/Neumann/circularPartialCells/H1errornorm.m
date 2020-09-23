%---------------------------------------------------------------------------%                          
%                           Error norm                                      %
%---------------------------------------------------------------------------%

function [err] = H1errornorm(ng,conn,elem,solution,flag,K,A,bc1,bc2,xC,yC,radius,flagsE);


%flag
%if 0 H0 norm
%if 1 H1 norm
%if 2 energy norm

[nelements, numberofnodes] = size(conn);

err  = 0;

for k =1:nelements
    
    if(flagsE(k,1) ~= 0)
    
    	el = elem(:,:,k);
    	%connEl = conn(k,:);
    	u_nE = elementSolution(solution,el);

    	for i =1:ng
        	for j=1:ng
            		[xi_g,eta_g,w1,w2] = gaussianquadrature(ng,i,j);
            		phi = shapefunction(xi_g, eta_g);
            		[x,y] = mapping(phi, el);
            
            		if(sqrt((xC- x)^2.0 + (yC - y)^2.0) <= radius)
                		
                                Dphi = shapefunctionFirstDerivative(xi_g,eta_g);
                		F = deformationGradient(el,Dphi);
                		jacob = jacobian(F);
                
                		if(flag == 0 || flag == 1)
                                        u = analyticalsolutionSinglePoint(K,A,bc1,bc2,x,y,xC,yC,radius);
                   			err = err + (phi*u_nE - u)^2*jacob*w1*w2; 
                		end
 
                		if(flag == 1)
                   			mat = u_nE'*Dphi*F^-1;
                   			[dudx,dudy] = derivativeAnalyticalsolution(K,A,x,y,xC,yC);
                   			err = err + ((mat(:,1) - dudx)^2 + (mat(:,2) - dudy)^2)*jacob*w1*w2; 
                		end

                                if(flag == 2)
                                        mat = u_nE'*Dphi*F^-1;
                                        [dudx,dudy] = derivativeAnalyticalsolution(K,A,x,y,xC,yC);
                                        err = err + A*((mat(:,1) - dudx)^2 + (mat(:,2) - dudy)^2)*jacob*w1*w2;
                                end

            		end 
        
         	end
    	end
    end
end

%take square root of everything
err = sqrt(abs(err));


%---------------------------------------------------------------------------%                          
%                           Error norm                                      %
%---------------------------------------------------------------------------%

function [err] = errornorms(ng,x_n,nn,solution,flag,K,A1,bc1,bc2,Lc,typeBC);


%flag
%if 0 H0 norm
%if 1 H1 norm

elem = elements(x_n,nn);
[nelements, numberofnodes] = size(elem);

err  = 0;


for k =1:nelements


        el = elem(k,:);
        u_nE = solution(k:(k+(nn-1)),1); 

        for i =1:ng

            [xi_g,w] = gaussianquadrature(ng,i);
            phi = shapefunction(xi_g, numberofnodes);
            x = mapping(phi, el);
            Dphi = shapefunctionFirstDerivative(xi_g,numberofnodes);
	    jacob = jacobian2(el,xi_g);
            doit = 1;
	    %If you are at the last element, and if the gauss point
	    %lies outside of the physical domain, do not use that
	    %point for integration. 
	    if((k == nelements) && (x> Lc))
		doit = 0;
	    end 
            if(doit == 1)
	    	if(flag == 0 || flag == 1)
	      		if (typeBC == 1)
				%analytical solution Dirichlet BC
				u = (-K/A1)*0.5 * x^2 + (((bc2 - bc1)/Lc) + (Lc/2)*(K/A1))*x + bc1;
	       		elseif(typeBC == 2)
				%analytical solution Neumann BC
				u = 0.5 * (-K/A1) * x^2 + (bc2 + (K/A1)* Lc) * x + bc1;
			end

      			err = err + (phi*u_nE - u)^2*(1/jacob)*w;
       		end

       		if(flag == 1)
			if (typeBC == 1)
				dudx = (-K/A1)*(x) + (((bc2 - bc1)/Lc) + (Lc/2)*(K/A1));
			elseif (typeBC == 2)
				%analytical solution Neumann BC
				dudx =  (-K/A1)*(x) + (bc2 + (K/A1)*Lc);
			end
       			mat  = u_nE'*Dphi'*jacob;
       			err = err + ((mat - dudx)^2)*(1/jacob)*w;
       		end
            end
       end
end


%take square root of everything
err = sqrt(abs(err));


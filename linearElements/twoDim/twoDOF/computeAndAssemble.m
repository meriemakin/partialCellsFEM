%---------------------------------------------------------------------------%                          
%                           Assembling of the element stiffness matrix      %
%                           and the element force vector                    % 
%---------------------------------------------------------------------------%


%elems are the elements of the cartesian grid
%flagsE: 0 if element outside of domain, 1 if element inside of domain 
%        2 if element partially inside of domain
%bcflags 1 if imposed to boundary conditions 0 if not
%dofconn is the degree of freedom connectivity
%tnn is the total number of nodes
%ng is the number of gauss points
%bc is the value of the boundary conditions on the boundary
function [s,f] = computeAndAssemble(elems,flagsE,bcflags,dofconn,tnn,ng,rhs,A,bc,xC, yC, radius);


[dim,nnpe,nel] = size(elems);
%s = zeros(tnn*2, tnn*2);
s = sparse(tnn*2, tnn*2);
f = zeros(tnn*2,1);

for i = 1:nel

    dofconnElem = dofconn(i,:);
    flag = flagsE(i,1);

    if(flag ~= 0)
    	sElem =  elementStiffness(elems(:,:,i),ng,nnpe,A,flag);
    	fElem = elementForceVector(ng,nnpe,elems(:,:,i),rhs,flag);
      
    	for j = 1:(nnpe*2)
                f(dofconnElem(j),1) = f(dofconnElem(j),1) + fElem(j,1);
        	for k = 1:(nnpe*2)
            		s(dofconnElem(j),dofconnElem(k)) = s(dofconnElem(j),dofconnElem(k)) + sElem(j,k);
        	end
    	end
   
    end
    
end

Ag = [];
bg = [];

for i=1:nel
    bcFlag = bcflags(i,1);
    [Al,bl] = boundaryInter(elems(:,:,i),bcFlag,bc,xC,yC,radius);
    dofconnElem = dofconn(i,:); 
    [count,dlt] = size(Al);
    for k=1:count
    	insertL = zeros(1,tnn*2);
    		for j = 1:(nnpe*2)
                     insertL(1,dofconnElem(j)) = Al(k,j);
    		end
    	insertR = bl(k,1);
    	Ag = [Ag;insertL];
    	bg = [bg;insertR];
    end
end


[count2,count3] = size(Ag);

for i=1:count2
    s(tnn*2-(i-1),:) = Ag(i,:);
    f(tnn*2-(i-1),1) = bg(i,1);
end

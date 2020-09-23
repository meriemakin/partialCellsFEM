function [x_n,y_n,elem] = neighborAveraging(flags, elem,conn,x_n,y_n,tnn,ncycles);


[r,c] = size(conn);

neighbornodes = cell(tnn-length(flags));
contrelements = cell(tnn-length(flags));
neighbor  = [];
contrelem = [];
pos = 1;

%go over all nodes
for i=1:tnn
    
    %treat nodes that do not lie on the domain boundary
    if(isempty(find(flags == i)))
       for j=1:r
         
         if(~isempty(find(conn(j,:) == i)))
            loc = find(conn(j,:) == i); 
            contrelem = [contrelem,j,loc];
            if(loc == 1)
               if(isempty(find(neighbor==conn(j,4))))
                 neighbor = [neighbor, conn(j,4)];
               end
               if(isempty(find(neighbor==conn(j,2))))
                 neighbor = [neighbor, conn(j,2)];
               end
            elseif(loc == 2)
               if(isempty(find(neighbor==conn(j,1))))
                 neighbor = [neighbor, conn(j,1)];
               end
               if(isempty(find(neighbor==conn(j,3))))
                 neighbor = [neighbor, conn(j,3)];
               end
            elseif(loc == 3)
               if(isempty(find(neighbor==conn(j,2))))
                 neighbor = [neighbor, conn(j,2)];
               end
               if(isempty(find(neighbor==conn(j,4))))
                 neighbor = [neighbor, conn(j,4)];
               end
            elseif(loc == 4)
               if(isempty(find(neighbor==conn(j,3))))
                 neighbor = [neighbor, conn(j,3)];
               end
               if(isempty(find(neighbor==conn(j,1))))
                 neighbor = [neighbor, conn(j,1)];
               end
            end
         end
      end
      neighbornodes(pos,1) = {[i, neighbor]};
      contrelements(pos,1) = {contrelem};
      neighbor = [];
      contrelem = [];
      pos = pos +1;
    end
end

%for i=1:pos-1
%disp('nodes')
%neighbornodes{i}
%disp('elems')
%contrelements{i}
%end

for i=1:ncycles
    for j=1:(pos-1)
         node = neighbornodes{j}(1);
         neighbors = neighbornodes{j}(2:end);

         contrelems = contrelements{j};
         x_n(node) = 0.0;
         y_n(node) = 0.0;
         
         for k=1:length(neighbors)
             x_n(node) = x_n(node) + x_n(neighbors(k));
             y_n(node) = y_n(node) + y_n(neighbors(k));
         end      
         x_n(node) = x_n(node)/(length(neighbors));
         y_n(node) = y_n(node)/(length(neighbors));
         
         for k=1:(0.5*length(contrelems))
             elem(1,contrelems(2*k),contrelems(2*k-1)) = x_n(node);
             elem(2,contrelems(2*k),contrelems(2*k-1)) = y_n(node);
         end         
    end

end
%figure
%plot(x_n,y_n,'+')


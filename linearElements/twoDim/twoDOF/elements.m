%---------------------------------------------------------------------------%                          
%                           generating elements from nodes                  %
%---------------------------------------------------------------------------%

%x_n and y_n are the nodes coordinates
%conn is the  connectivity 

function elem = elements(x_n, y_n,conn);

[n,l] = size(conn);
elem = zeros(2,l,n);

for i =1:n
    connElem = conn(i,:);
    elem(1,1:l,i) = [x_n(connElem(1:l))];
    elem(2,1:l,i) = [y_n(connElem(1:l))];
end

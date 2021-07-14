function k = elemStiff(elemCoord, materials, problem)
%{
This function takes the coordinates of the element (elemCoord matrix), a
vector of the element material properties [E, nu, t], and the problem type
('plane stress / plane strain) as the input arguments and the conduct the
numerical integration to calculate of the element stiffness matrix
%}

[nElemNode,nDim] = size(elemCoord);

% extracting the associated material properties
E  = materials(1);
nu = materials(2);
t  = materials(3);

% calling the "strainStressMatrix" function to calculate the [D] matrix
D = strainStressMatrix (E, nu, problem);


% generating a list of gauss quadrature integration points and weights
% (in natural coordinates)
[xiList, wList] = integrationPoints(nDim,nElemNode);


% initializing the variable to store the [B] and [k] matrices values
B = zeros(3, nElemNode*2);
k = zeros(2*nElemNode);

% loop over all integration points to conduct the numerical integration and
% obtaining the element stiffness matrix
switch nDim
    % nPoints = number of the integration points
    %=====================================================================%
    case 1 % 1D
        %-----------------------------------------------------------------%    
        if nElemNode == 2
           J = (elemCoord(2)-elemCoord(1))/2;
           [N, dNdxi] = shapeFunctions(nDim,nElemNode,xiList);
           B(1,1) = dNdxi(1)/J;
           B(1,3) = dNdxi(2)/J;
           k = B'*D*B*J*t*wList(1) + B'*D*B*J*t*wList(2);
        %-----------------------------------------------------------------%    
        elseif nElemNode == 3
           nPoints = 3;
           for n=1:nPoints
           B = zeros(3, nElemNode*2);
           [N, dNdxi] = shapeFunctions(nDim,nElemNode,xiList(n))
               J= (elemCoord(1) - 2*elemCoord(2) + elemCoord(3))*xiList(n) - elemCoord(1)/2 + elemCoord(3)/2;
               for i=1:nElemNode
               B(1,2*i-1) = dNdxi(i)/J;
               end
           k = k+ B'*D*B*J*t*wList(n);
           end
        
         
        %-----------------------------------------------------------------%    
        end
    %=====================================================================%    
    case 2 % 2D
        %-----------------------------------------------------------------%
        if (nElemNode == 2) || (nElemNode == 3)
            [N, dNdxi] = shapeFunctions(nDim,nElemNode,xiList(:,1));
            J = [elemCoord(1,1)-elemCoord(3,1) elemCoord(1,2)-elemCoord(3,2);
                elemCoord(2,1)-elemCoord(3,1) elemCoord(2,2)-elemCoord(3,2)];
            dNdx = J^-1*dNdxi';
            for i=1:nElemNode
                B(1,2*i-1) = dNdx(1,i);
                B(2,2*i) = dNdx(2,i);
                B(3,2*i-1)= dNdx(2,i); B(3,2*i)= dNdx(1,i);
            end
          k =  B'*D*B*det(J)*t*wList(1); 
                
            
        %-----------------------------------------------------------------%    
        elseif nElemNode == 6
            nPoints = 3;
            for i = 1:nPoints
            for n=1:nPoints
                list =[];
                list(1)=xiList(1,i); list(2)=xiList(2,n);
            [N, dNdxi] = shapeFunctions(nDim,nElemNode,list);
            J = [elemCoord(3,1)*(2*xiList(1,i)+2*xiList(2,n)-1)+2*elemCoord(3,1)*(xiList(1,i)+xiList(2,n)-1)-4*elemCoord(6,1)*(xiList(1,i)+xiList(2,n)-1)+elemCoord(1,1)*(4*xiList(1,i)-1)-4*xiList(1,i)*elemCoord(6,1)+4*xiList(2,n)*elemCoord(4,1)-4*xiList(2,n)*elemCoord(5,1)  
                 elemCoord(3,2)*(2*xiList(1,i)+2*xiList(2,n)-1)+2*elemCoord(3,2)*(xiList(1,i)+xiList(2,n)-1)-4*elemCoord(6,2)*(xiList(1,i)+xiList(2,n)-1)+elemCoord(1,2)*(4*xiList(1,i)-1)-4*xiList(1,i)*elemCoord(6,2)+4*xiList(2,n)*elemCoord(4,2)-4*xiList(2,n)*elemCoord(5,2);
                 elemCoord(3,1)*(2*xiList(1,i)+2*xiList(2,n)-1)-elemCoord(5,1)*(4*xiList(1,i)+4*xiList(2,n)-4)+2*elemCoord(3,1)*(xiList(1,i)+xiList(2,n)-1)+elemCoord(2,1)*(4*xiList(2,i)-1)+4*xiList(1,i)*elemCoord(4,1)-4*xiList(1,i)*elemCoord(6,1)-4*xiList(2,n)*elemCoord(5,1) 
                 elemCoord(3,2)*(2*xiList(1,i)+2*xiList(2,n)-1)-elemCoord(5,1)*(4*xiList(1,i)+4*xiList(2,n)-4)+2*elemCoord(3,2)*(xiList(1,i)+xiList(2,n)-1)+elemCoord(2,2)*(4*xiList(2,i)-1)+4*xiList(1,i)*elemCoord(4,2)-4*xiList(1,i)*elemCoord(6,2)-4*xiList(2,n)*elemCoord(5,2)]
            dNdx = J^-1*dNdxi';
            for j=1:nElemNode
                B(1,2*j-1) = dNdx(1,j);
                B(2,2*j) = dNdx(2,j);
                B(3,2*j-1)= dNdx(2,j); B(3,2*j)= dNdx(1,j);
            end
            k = k+B'*D*B*det(J)*t*wList(n)*wList(i);
            end
            end
            
        %-----------------------------------------------------------------%    
        elseif nElemNode == 4
            nPoints = 4;
            for i = 1:nPoints
            for n = 1:nPoints
                list =[];
                list(1)=xiList(1,i); list(2)=xiList(2,n);
            [N, dNdxi] = shapeFunctions(nDim,nElemNode,list);
            
            J = 1/4*[-(1-xiList(2,n))*elemCoord(1,1)+(1-xiList(2,n))*elemCoord(2,1)+(1+xiList(2,n))*elemCoord(3,1)-(1+xiList(2,n))*elemCoord(4,1)  -(1-xiList(2,n))*elemCoord(1,2)+(1-xiList(2,n))*elemCoord(2,2)+(1+xiList(2,n))*elemCoord(3,2)-(1+xiList(2,n))*elemCoord(4,2);
                        -(1-xiList(1,i))*elemCoord(1,1)-(1+xiList(1,i))*elemCoord(2,1)+(1+xiList(1,i))*elemCoord(3,1)+(1-xiList(1,i))*elemCoord(4,1)  -(1-xiList(1,i))*elemCoord(1,2)-(1+xiList(1,i))*elemCoord(2,2)+(1+xiList(1,i))*elemCoord(3,2)+(1-xiList(1,i))*elemCoord(4,2)];
            dNdx = J^-1*dNdxi';
            for j=1:nElemNode
                B(1,2*j-1) = dNdx(1,j);
                B(2,2*j) = dNdx(2,j);
                B(3,2*j-1)= dNdx(2,j); B(3,2*j)= dNdx(1,j);
            end
            
            k = k + B'*D*B*det(J)*t*wList(n)*wList(i);  
            end
            end
        %-----------------------------------------------------------------%    
        elseif nElemNode == 8
            nPoints = 9;
            for i=1:nPoints
            for n=1:nPoints
                list =[];
                list(1)=xiList(1,i); list(2)=xiList(2,n);
            [N, dNdxi] = shapeFunctions(nDim,nElemNode,list);
            J = [elemCoord(8,1)*(xiList(2,n)^2-1)/2-elemCoord(6,1)*(xiList(2,n)^2-1)/2-elemCoord(1,1)*(xiList(1,i)-1)*(xiList(2,n)-1)/4-elemCoord(2,1)*(xiList(1,i)+1)*(xiList(2,n)-1)/4+elemCoord(3,1)*(xiList(1,i)+1)*(xiList(2,n)+1)/4+elemCoord(4,1)*(xiList(1,i)-1)*(xiList(2,n)+1)/4-elemCoord(1,1)*(xiList(2,n)-1)*(xiList(1,i)+xiList(2,n)+1)/4+elemCoord(3,1)*(xiList(2,n)+1)*(xiList(1,i)+xiList(2,n)-1)/4+xiList(1,i)*elemCoord(5,1)*(xiList(2,n)-1)-xiList(1,i)*elemCoord(7,1)*(xiList(2,n)+1)+elemCoord(2,1)*(xiList(2,n)-1)*(xiList(2,n)-xiList(1,i)+1)/4+elemCoord(4,1)*(xiList(2,n)+1)*(xiList(1,i)-xiList(2,n)+1)/4
                elemCoord(8,2)*(xiList(2,n)^2-1)/2-elemCoord(6,2)*(xiList(2,n)^2-1)/2-elemCoord(1,2)*(xiList(1,i)-1)*(xiList(2,n)-1)/4-elemCoord(2,2)*(xiList(1,i)+1)*(xiList(2,n)-1)/4+elemCoord(3,2)*(xiList(1,i)+1)*(xiList(2,n)+1)/4+elemCoord(4,2)*(xiList(1,i)-1)*(xiList(2,n)+1)/4-elemCoord(1,2)*(xiList(2,n)-1)*(xiList(1,i)+xiList(2,n)+1)/4+elemCoord(3,2)*(xiList(2,n)+1)*(xiList(1,i)+xiList(2,n)-1)/4+xiList(1,i)*elemCoord(5,2)*(xiList(2,n)-1)-xiList(1,i)*elemCoord(7,2)*(xiList(2,n)+1)+elemCoord(2,2)*(xiList(2,n)-1)*(xiList(2,n)-xiList(1,i)+1)/4+elemCoord(4,2)*(xiList(2,n)+1)*(xiList(1,i)-xiList(2,n)+1)/4;
                elemCoord(5,1)*(xiList(1,i)^2-1)/2-elemCoord(7,1)*(xiList(1,i)^2-1)/2-elemCoord(1,1)*(xiList(1,i)-1)*(xiList(2,n)-1)/4+elemCoord(2,1)*(xiList(1,i)+1)*(xiList(2,n)-1)/4+elemCoord(3,1)*(xiList(1,i)+1)*(xiList(2,n)+1)/4-elemCoord(4,1)*(xiList(1,i)-1)*(xiList(2,n)+1)/4-elemCoord(1,1)*(xiList(1,i)-1)*(xiList(1,i)+xiList(2,n)+1)/4+elemCoord(3,1)*(xiList(1,i)+1)*(xiList(1,i)+xiList(2,n)-1)/4-xiList(2,n)*elemCoord(6,1)*(xiList(1,i)+1)+xiList(2,n)*elemCoord(8,1)*(xiList(1,i)-1)+elemCoord(2,1)*(xiList(1,i)+1)*(xiList(2,n)-xiList(1,i)+1)/4+elemCoord(4,1)*(xiList(1,i)-1)*(xiList(1,i)-xiList(2,n)+1)/4
                elemCoord(5,2)*(xiList(1,i)^2-1)/2-elemCoord(7,2)*(xiList(1,i)^2-1)/2-elemCoord(1,2)*(xiList(1,i)-1)*(xiList(2,n)-1)/4+elemCoord(2,2)*(xiList(1,i)+1)*(xiList(2,n)-1)/4+elemCoord(3,2)*(xiList(1,i)+1)*(xiList(2,n)+1)/4-elemCoord(4,2)*(xiList(1,i)-1)*(xiList(2,n)+1)/4-elemCoord(1,2)*(xiList(1,i)-1)*(xiList(1,i)+xiList(2,n)+1)/4+elemCoord(3,2)*(xiList(1,i)+1)*(xiList(1,i)+xiList(2,n)-1)/4-xiList(2,n)*elemCoord(6,2)*(xiList(1,i)+1)+xiList(2,n)*elemCoord(8,2)*(xiList(1,i)-1)+elemCoord(2,2)*(xiList(1,i)+1)*(xiList(2,n)-xiList(1,i)+1)/4+elemCoord(4,2)*(xiList(1,i)-1)*(xiList(1,i)-xiList(2,n)+1)/4]
            dNdx = J^-1*dNdxi';
            for j=1:nElemNode
                B(1,2*j-1) = dNdx(1,j);
                B(2,2*j) = dNdx(2,j);
                B(3,2*j-1)= dNdx(2,j); B(3,2*j)= dNdx(1,j);
            end
            k = k+B'*D*B*det(J)*t*wList(n)*wList(i); 
            end
            end
            
        %-----------------------------------------------------------------%   
        end   
    %=====================================================================%
end




















end

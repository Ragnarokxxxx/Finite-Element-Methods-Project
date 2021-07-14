function r = elemEquivalentLoad(edgeCoord, load)
%{
This function takes the coordinates of an edge of an element, and the
componenets of the applied traction force [Tx, Ty], and calculates the
equivalent nodal forces by conducting numerical integration
%}

[nEdgeNode, nDim] = size(edgeCoord);
nDof = length(load);

% initializing the variable to store the equivalent nodal forces
r = zeros(nDof*nEdgeNode, 1);

% calling the "integrationPoints" function to generate a list of
% Guass quadrature integration points and weights
[xiList, wList] = integrationPoints(nDim-1,nEdgeNode);



% loop over integration points to conduct the numerical integration and
% calculate the equivalent nodal forces

        %-----------------------------------------------------------------%    
        if nEdgeNode == 2
            nPoints = 2;
            Tx = load(1); Ty = load(2);
            T = [Tx;Ty];
            f = zeros(2);
            
            for n=1:nPoints
                [N, dNdxi] = shapeFunctions(nDim-1,nEdgeNode,xiList(n));
                J = ((edgeCoord(2,2)-edgeCoord(1,2))^2+(edgeCoord(2,1)-edgeCoord(1,1))^2)/2;
                f = f+T*N'*J*wList(n);
            end    
            for n=1:nEdgeNode
                r(2*n-1,1) =  f(1,n);
                r(2*n,1)= f(2,n);
            end
        %-----------------------------------------------------------------%    
        elseif nEdgeNode == 3
            nPoints = 3;
            Tx = load(1); Ty = load(2);
            T = [Tx;Ty];
            f = zeros(2,3);
            for n=1:nPoints
                [N, dNdxi] = shapeFunctions(nDim-1,nEdgeNode,xiList(n));
                J = ((edgeCoord(3,2)-edgeCoord(1,2))^2+(edgeCoord(3,1)-edgeCoord(1,1))^2)/2;
                f = f+T*N'*J*wList(n);
            end    
            for n=1:nEdgeNode
                r(2*n-1,1) =  f(1,n);
                r(2*n,1)= f(2,n);
            end
            
        %-----------------------------------------------------------------%    
        end
    %=====================================================================%    
    
    %=====================================================================%

    %=====================================================================%















end


function [xi, w] = integrationPoints(nDim,nElemNode)
%{
This function takes the nDim and nElemNode as inputs and returns the
following outputs:

(1) xi: a (nDim-by-nIntegrationPoints) matrix including the gauss
quadrature points

(2) w:  a list (vector) of gauss quadrature weights
%}

switch nDim
    % nPoints = number of the integration points
    %=====================================================================%
    case 1 % 1D
        nPoints = nElemNode;
        %-----------------------------------------------------------------%
        if nElemNode == 1
            xi = 0.0;
            w  = 2.0;
        %-----------------------------------------------------------------%    
        elseif nElemNode == 2
            xi = [-1/sqrt(3), 1/sqrt(3)];
            w  = [1.0       , 1.0      ];
        %-----------------------------------------------------------------%    
        elseif nElemNode == 3
            % TODO
            % xi = ...
            xi = [-sqrt(3/5) 0 sqrt(3/5)];
            % w  = ...
            w = [5/9 8/9 5/9];
            
        %-----------------------------------------------------------------%    
        end
    %=====================================================================%    
    case 2 % 2D
        %-----------------------------------------------------------------%
        if (nElemNode == 2) || (nElemNode == 3)
            nPoints = 1;
            
            % TODO
            % xi = ...
            xi = [1/3; 1/3];
            % w  = ...
            w = [1/2];
            
        %-----------------------------------------------------------------%    
        elseif nElemNode == 6
            nPoints = 3;
            
            % TODO
            % xi = ...
            xi = [1/2 0 1/2; 1/2 1/2 0];
            % w  = ...
            w = [1/6 1/6 1/6];
            
        %-----------------------------------------------------------------%    
        elseif nElemNode == 4
            nPoints = 4;
            
            % TODO
            % xi = ...
            xi = [-sqrt(1/3) sqrt(1/3) -sqrt(1/3) sqrt(1/3);
                  -sqrt(1/3) -sqrt(1/3) sqrt(1/3) sqrt(1/3)];
            % w  = ...
            w = [1 1 1 1];
        %-----------------------------------------------------------------%    
        elseif nElemNode == 8
            nPoints = 9;
            
            % TODO
            % xi = ...
            xi = [-sqrt(3/5) 0 sqrt(3/5) -sqrt(3/5) 0
                  sqrt(3/5) -sqrt(3/5) 0 sqrt(3/5);
                  -sqrt(3/5) -sqrt(3/5) -sqrt(3/5)
                  0 0 0 sqrt(3/5) sqrt(3/5) sqrt(3/5)];
            % w  = ...
            w = [25/81 40/81 25/81 40/81 64/81 40/81 25/81 40/81 25/81];
            
        %-----------------------------------------------------------------%   
        end   
    %=====================================================================%
end

end


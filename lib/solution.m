function [displacements, reactionForces] = solution(K, R, bcIndex, bcValues)
%{
This function imposes the boundary conditions on the assembled global
matrices and calculates the displacements and reaction forces
%}

dofTotal = length(R);                     % number of all dofs
dofActive = setdiff(1:dofTotal, bcIndex); % a list (vector) of free dofs



% imposing boundary conditions
% TODO
% Kf = ...
Kf = zeros(length(dofActive));
i_Kf = 1;
for i_Dof = 1:length(K)
    j_Kf = 1;
    if any(dofActive == i_Dof)
        for j_Dof = 1:length(K)
            if any(dofActive == j_Dof)
                abc = K(i_Dof,j_Dof);
                Kf(i_Kf,j_Kf) = Kf(i_Kf,j_Kf)+abc;
                j_Kf = j_Kf + 1;
            end
        end
        i_Kf = i_Kf + 1;
    end
end

        
% Rf = ...
Rf = R(dofActive);

% solving the system
% TODO
% u =...
u = Kf^-1*Rf; 


% forming the displacements vector using both free dofs values and the
% prescribed displacements values

% TODO
displacements            = zeros(dofTotal,1);
% displacements = ...
displacements(dofActive) = u;



% calculating reaction forces

% TODO
% reactionForces = ...
reactionForces = K*displacements-R;
end


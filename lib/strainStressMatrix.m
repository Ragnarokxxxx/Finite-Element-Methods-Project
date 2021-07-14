function D = strainStressMatrix (E, nu, problem) 
%{
This function takes material properties (E and nu), and the problem type
(plain stress/ plane strain) as the inputs and generates the strain-stress
matrix ([D])
%}

problem = upper(erase(problem,' '));


switch problem
    case 'PLANESTRESS'
        % TODO
        % D = ...
        D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    
    case 'PLANESTRAIN'
        % TODO
        % D = ...
        D = E/((1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-nu)/2];
end

end

function printOutput(Results)
    displacements  = Results.displacements;
    reactionForces = Results.reactions;
    stresses       = Results.stresses;
    strains        = Results.strains;
    naturalCoord   = Results.naturalcoordinates;
    globalCoord    = Results.globalcoordinates;

    % TODO
    % display the results on command window
    fprintf('==================================================================\n')
    fprintf('                   O U T P U T    S U M M A R Y                   \n')
    fprintf('==================================================================\n')
    fprintf('------------------------------------------------------------------\n')
    fprintf('                        Nodal Displacements                       \n')
    fprintf('------------------------------------------------------------------\n')
    for i = 1:size(displacements,1)
        fprintf('[%3i]%14.5e%14.5e\n', i,displacements(i,1),displacements(i,2))
    end
    fprintf('------------------------------------------------------------------\n')
    fprintf('                         Reaction Forces                          \n')
    fprintf('------------------------------------------------------------------\n')
    for i = 1:size(reactionForces,1)
        fprintf('[%3i]%15.5f\n', i,reactionForces(i))
    end
    fprintf('-------------------------------------------------------------------------------------------\n')
    fprintf('                                    Strains and Stresses                                   \n')
    fprintf('-------------------------------------------------------------------------------------------\n')
    for iElem = 1:size(stresses,3)
        fprintf(['                             Element Number = ' num2str(iElem) '                           \n'])
        fprintf('----------------------------------------------------------------------------------------------------------------------------------------------\n')
        fprintf('Intergration\t\tLocation\t\t\t\t\t\t\tStrains\t\t\t\t\t\t\t\t\tStresses\n')
        fprintf('\tPoint\t\t\t(Natural)\t\t\t\t\te11\t\t\te22\t\t\te12\t\t\t\ts11\t\t\ts22\t\t\ts12\n')
        fprintf('----------------------------------------------------------------------------------------------------------------------------------------------\n')
        for iPoints = 1:size(naturalCoord,2)
            fprintf('%i\t\t\t\t(%8.4f%8.4f)\t\t(%8.4f%8.4f)\t\t', iPoints, [naturalCoord(1,iPoints,iElem),naturalCoord(2,iPoints,iElem)])
            fprintf('%13.4e%13.4e%13.4e\t', strains(1,iPoints,iElem), strains(2,iPoints,iElem), strains(3,iPoints,iElem))
            fprintf('%13.4e%13.4e%13.4e\t\n', stresses(1,iPoints,iElem), stresses(2,iPoints,iElem), stresses(3,iPoints,iElem))
        end
        fprintf('\n\n')
    end
    fprintf('=============================================================================================================================================\n')
    
end



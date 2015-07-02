function [pairsWithinGroups, pairsBetweenGroups] = getPairsWithinAndBetweenGroups(groupIds, okPairs, maxPairs)
    %input: list of group ids: eg [10 10 10, 12, 12, 13, 14, 14, 14, 14]
    %output: pairsWithinGroups: [1 2; 1 3; 2 3; 4 5; etc]
           % pairsBetweenGroups: [1 4; 1 5; etc]

    useMax = false;
           
    if ~exist('maxPairs', 'var') || isempty(maxPairs)
        maxPairs = 1e8;
    end
    
    
    [uniqueGroups, groupCounts] = uniqueCount(groupIds);
    nGroups = length(uniqueGroups);
    % renumber groupIds from 1 : nGroups    
    
    groupIds_renum = zeros(size(groupIds));
    for i = 1:nGroups
        groupIds_renum(groupIds==uniqueGroups(i)) = i;
    end
    groupIds = groupIds_renum;
    
    % find number of within-group pairs
    nPairsWithinGroups = sum ( (groupCounts .* (groupCounts-1)) / 2);
    if useMax && (nPairsWithinGroups > maxPairs)
        fprintf(' [Found %d pairs - capping at max of %d]\n', nPairsWithinGroups, maxPairs);        
        nPairsWithinGroups = maxPairs;
    end
    pairsWithinGroups  = zeros(nPairsWithinGroups,2);
    pw_i = 1;

    % find number of between-group pairs
    M = groupCounts(:) * groupCounts(:)';
    for i = 1:length(groupCounts)
        M(i,i) = 0;
    end
	nPairsBetweenGroups = sum(M(:))/2;
    if useMax
        nPairsBetweenGroups = min(nPairsBetweenGroups, maxPairs);
    end
    pairsBetweenGroups  = zeros(nPairsBetweenGroups,2);
    pb_i = 1;
    
    Ntotal = nPairsBetweenGroups + nPairsBetweenGroups;
    
    % check if validating between-group pairs.
    validatePairs = exist('okPairs', 'var') && ~isempty(okPairs);

    for i = 1:length(groupIds)        
        for j = i+1:length(groupIds)
%             progressBar;
            
            GrpI = groupIds(i);
            GrpJ = groupIds(j);
                        
            if (GrpI == GrpJ) 
                if (~useMax) || (pw_i <= maxPairs)
                    if pw_i > size(pairsWithinGroups,1)
                        3;
                    end
                    pairsWithinGroups(pw_i,1:2) = [i j];

                    pw_i = pw_i + 1;
                end
                
            else % => (GrpI ~= GrpJ) 
                if ((~useMax) || (pb_i <= maxPairs)) && (~validatePairs || okPairs(GrpI, GrpJ))
                    if pb_i > size(pairsBetweenGroups,1)
                        3;
                    end
                    pairsBetweenGroups(pb_i,1:2) = [i j]; 
                    pb_i = pb_i + 1;                
                end
            end
            
        end        
        if (pw_i > maxPairs) && (pb_i > maxPairs)
            break;
        end
    end
%     progressBar('done');
    
    pairsWithinGroups (pw_i:end,:) = [];
    pairsBetweenGroups(pb_i:end,:) = [];
end
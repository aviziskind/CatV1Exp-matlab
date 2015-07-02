function [N, vals, idx] = binCountForPairs_Matlab(pairIdxs, tf, binIds, catIds, nBins, nCats, allVals)
    getVals = nargout > 1;
    getIdx  = nargout > 2;
    
    if getVals
        vals = zeros(1,length(pairIdxs));
        count = 1;
    end
    if getIdx
        idx = zeros(1,length(pairIdxs));
    end
    

    N = zeros(nBins, nCats);
    if ~isempty(catIds)        
        for i = pairIdxs(:)'
            if tf(i) && binIds(i) && catIds(i)  % ie. binId & catId > 0
                N(binIds(i), catIds(i)) = N(binIds(i), catIds(i))+1;
                if getIdx
                    idx(count) = i;
                end
                if getVals
                    vals(count) = allVals(i);
                    count = count +1;
                end
            end
        end    
    else            
        for i = pairIdxs(:)'
            if tf(i) && binIds(i)  % ie. binId  > 0
                N(binIds(i)) = N(binIds(i))+1;
                if getIdx
                    idx(count) = i;
                end                
                if getVals
                    vals(count) = allVals(i);
                    count = count +1;
                end
            end
        end    
    end

    if getVals
        vals = vals(1:count-1);
    end
    if getIdx
        idx = idx(1:count-1);
    end
    
end
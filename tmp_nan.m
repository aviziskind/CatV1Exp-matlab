fn = fieldnames(allSpkStimHists.(psthType));

for i = 1:length(fn)
    
    psth_name = fn{i};
    
    s = allSpkStimHists.(psthType).(psth_name);        
    [bins, uVals, meanRate, bckgRate] = deal(...
        s.bins, s.uVals, s.meanRate, s.bckgRate);
    if any(isnan(uVals))
    
        if isfield(s, 'uValsLists')
            uValLists = s.uValsLists;
            allHistVals = zeros(s.dims);
            for i = 1:length(uVals)
                allHistVals(uValLists{i}) = uVals(i);
            end            
        elseif isfield(s, 'vals_uint8')
            allHistVals_uint8 = s.vals_uint8;
            allHistVals = uVals(allHistVals_uint8);            
        end
        
        
        [uVals, uValsLists] = uniqueList(allHistVals(:));   
        valCounts = cellfun(@length, uValsLists);
        assert(uVals(1) == 0);            
        n_total = numel(allHistVals);            
        n_rest = sum( valCounts(2:end) );
        if (n_rest*4) < n_total  % space for listing indices of non zero entries (using 4-byte integers) < space for entries as 1-byte value indices.
            uValsLists = cellfun(@uint32, uValsLists, 'un', 0);                
            allSpkStimHists.(psthType).(psth_name) = struct('bins', bins, 'uVals', uVals(2:end), 'uValsLists', {uValsLists(2:end)}, 'dims', size(allHistVals), 'meanRate', meanRate, 'bckgRate', bckgRate);
        else                   % mostly non-zeros: store indices to data values                
            allHistVals_uint8 = zeros(size(allHistVals), 'uint8');
            for val_i = 1:length(uVals)
                allHistVals_uint8(uValsLists{val_i}) = val_i;
            end            
            allSpkStimHists.(psthType).(psth_name) = struct('bins', bins, 'uVals', uVals, 'vals_uint8', allHistVals_uint8, 'meanRate', meanRate, 'bckgRate', bckgRate);
        end                        
        
        
        
        
        
    end

    
end
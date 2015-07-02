function dbFixTempFreqBatchStimIds
    % there were some "temporal frequency batch" drifting grating experiments -
    % experiments where temporal frequency was varied in different
    % presentations, but their stimulus id (supposed to be 7) is
    % mislabeled (all as "9"). This script fixes that.
    
    hnd = dbOpenExpDb;

    % Find all drifting grating experiments where are multiple temporal
    % frequencies.

%   {'Single Grating', 'Flashed Grating Batch', 'Orientation Batch', 'Spatial Frequency Batch', 'Temporal Frequency Batch', 'Free Grating Batch'}
%    1, 3, 2, 6, 7, 9, 10

    allDids = unique(getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_GRATING_PRES'));
    nTempFreqs = zeros(1, length(allDids));    
    isTempFreqBatch = false(1, length(allDids));    
    
    progressBar('init-', length(allDids));
    for i = 1:length(allDids)
        progressBar(i);
        tmpFreqs = getFieldsFromDatabaseTable(hnd, 'DBL_TEMP_PERIOD_FRM', 'TBL_GRATING_PRES', {'DATAFILE_ID', allDids(i)});
        nTempFreqs(i) = length(unique(tmpFreqs));
        if nTempFreqs > 1
            isTempFreqBatch(i) = true;
        end       
    end
    progressBar('done');
    
    allStimIds = getFieldsFromDatabaseTable(hnd, 'STIMULUS_TYPE_ID', 'TBL_DATA_FILES');
    
    nTFbatch = nnz(isTempFreqBatch);
    tfBatchInds = find(isTempFreqBatch);
    tfbStimIds = zeros(1,nTFbatch);    
    for i = 1:nTFbatch
        did = allDids(tfBatchInds(i));
        tfbStimIds(i) = getFieldsFromDatabaseTable(hnd, 'STIMULUS_TYPE_ID', 'TBL_DATA_FILES', {'DATAFILE_ID', did});
    end
    
    nIncorrect = nnz(tfbStimIds(i) ~= 7);        
    fprintf('There are %d TEMPORAL FREQUENCY BATCH grating experiments, and %d are mislabeled.\n', nTFbatch, nIncorrect);
    if nIncorrect > 0    
        nCorrected = 0;
        for i = 1:nTFbatch
            if (tfbStimIds(i) ~= 7)
                did = allDids(tfBatchInds(i));
                updateValueInDatabaseTable(hnd, 7, 'STIMULUS_TYPE_ID', 'TBL_DATA_FILES', {'DATAFILE_ID', did});
                nCorrected = nCorrected + 1;
            end        
        end
        fprintf('Corrected %d entries.', nCorrected);
    end
    
    
end


% varBreakdown(nTempFreqs(isTempFreqBatch))
%     # Exp     Number of _temp_freqs
%  ---------   -----------------------
%       18   :    2 
%       2    :    3 
%       2    :    4 
%       1    :    5 
%       2    :    6 
%       9    :    7 
%       4    :    8 


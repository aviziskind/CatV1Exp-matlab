% function generateAllNoiseStims

    hnd = dbOpenExpDb;
    tableName = 'TBL_NOISE_PRES';
    allNGradations = unique(getFieldsFromDatabaseTable(hnd, 'LNG_N_GRADATIONS', tableName))';
    
    for nGrad = allNGradations % [2 3 8 256]
%         criteria = {'LNG_N_MISSED_FRM', 0; 'LNG_N_GRADATIONS', nGrad};
        criteria = {'LNG_N_GRADATIONS', nGrad};
        [nSustained, nDisplayed, ncols, nrows, randSeed] = getFieldsFromDatabaseTable(hnd, ... 
            {'LNG_N_SUSTAINED_FRM', 'LNG_N_DISPLAYED_FRM', 'LNG_N_COLUMNS', 'LNG_N_ROWS', 'LNG_RAND_SEED'}, tableName, criteria);
        
        nFrames = max([nSustained nDisplayed],[],2);
%         nPixels1 = ncols .* nrows .* nSustained;
%         nPixels2 = ncols .* nrows .* nDisplayed;
        [nMaxPixels ind] = max(  nrows.*ncols.* nFrames  );
        
        disp(['Generating noise movie file for nGrad = ' num2str(nGrad) '.']);
        generateNoiseStim([nrows(ind) ncols(ind)], nFrames(ind), nGrad, randSeed(ind)) ;
    end
        
        
%         criteria = {'LNG_N_MISSED_FRM', 0; 'LNG_N_GRADATIONS', nGrad};

%     disp(['Generating noise movie file for nGrad = ' num2str(nGrad) '.']);
%     generateNoiseStim(Did(ind), noisePresId(ind), 'saveRaw');
    
%     sortrows([nPixels1 / 1e6, nPixels2 / 1e6, Did, nGradations])

    
    
%     % shows that don't save any disk space by using nrow_ncol instead of
%     % npix:
% 
%     criteria = {'LNG_N_MISSED_FRM', 0};
%     [Did, noisePresId, nSustained, nDisplayed, ncols, nrows] = getFieldsFromDatabaseTable(hnd, ... 
%         {'DATAFILE_ID', 'NOISE_PRES_ID', 'LNG_N_SUSTAINED_FRM', 'LNG_N_DISPLAYED_FRM', 'LNG_N_COLUMNS', 'LNG_N_ROWS'}, tableName, criteria);
% 
%     npixels = ncols .* nrows;
%     [b,m,n] = unique([ncols nrows], 'rows');
%     A = [ncols(m), nrows(m), npixels(m)]
%     ind = ord( npixels(m) );
%     Data = [npixels(m), ncols(m), nrows(m)];
%     Data = Data(ind, :)

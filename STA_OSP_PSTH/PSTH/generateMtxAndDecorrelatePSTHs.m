function r = generateMtxAndDecorrelatePSTHs(R, frameStimIds, nBinsPerFrame, nFramesPerExtFrame, r_framesToConsider)

    % 1. Generates the matrix for PSTH decorrelation
        
    nTotalFrames = length(frameStimIds);
    stimIds = setdiff( unique(frameStimIds), 0); % can have some 0 values, if are ignoring certain stimuli.
    nStims = length( stimIds );
    assert( all ( stimIds(:)' == 1:nStims ));
    if ~exist('r_framesToConsider', 'var') || isempty(r_framesToConsider)    
        r_framesToConsider = 1:nFramesPerExtFrame;        
    end
    R_framesToConsider = 1:nFramesPerExtFrame; % consider all R values (-> rectangular matrix if not considering all r values)
%     R_framesToConsider = r_framesToConsider; % consider only R values that correspond to considered r values
    
    nrFramesConsidered = length(r_framesToConsider);
    nRFramesConsidered = length(R_framesToConsider);
    nR = nStims * nRFramesConsidered;
    nr = nStims * nrFramesConsidered;
    matDims = [nR, nr];
    matSize = nR*nr;
    
    printDetails = true;
    usePrecision = 'single';
    
    if isempty(usePrecision)        
        % check memory allocation.
        userview = memory;
        maxArrayBytes = userview.MaxPossibleArrayBytes; % my machine: 927580160 (885 MB)
        nMaxDouble = floor(maxArrayBytes/8)-1;
        nMaxSingle = floor(maxArrayBytes/4)-1;
        nMaxInt8  =  floor(maxArrayBytes)  -1;
        if matSize < nMaxDouble       % on my machine: enough for 10766 ^2 mtx
            usePrecision = 'double';
        elseif matSize < nMaxSingle   % on my machine: enough for 15227 ^2 mtx
            usePrecision = 'single';
        elseif matSize < nMaxInt8     % on my machine: enough for 30455 ^2 mtx
            usePrecision = 'uint8';
        else
            error('Matrix is too big. We can start off with a sparse matrix, but this will be very slow');
        end
    end
    
    A = zeros(matDims, usePrecision);
    if printDetails
        d = whos('A');
        siz_mb = d.bytes / (1024^2);
        fprintf(['Using ' usePrecision ' precision. Matrix is %d x %d (%.0f MB)'], nR, nr, siz_mb);
    end
    
    showWorking = false + 0;    
    if showWorking 
        showDetails = all(matDims < 50);
        figure(162);
        hM = imagesc(A);
        set(gca, 'tickdir', 'out');
        caxis([0  length(find(frameStimIds == 1))]);
        colormap('jet');
        
        if showDetails        
            drawHorizontalLine(.5 + [1:nR-1], 'Linestyle', ':');
            drawVerticalLine(  .5 + [1:nr-1], 'Linestyle', ':');
            drawHorizontalLine(.5 + [1:nStims] * nRFramesConsidered, 'LineWidth', 1 );
            drawVerticalLine(  .5 + [1:nStims] * nrFramesConsidered, 'LineWidth', 1 );
            if nr == nR
                axis square;
            end

            stringsh = zeros(matDims);
            for i = 1:nR
                for j = 1:nr
                    stringsh(i,j) = text(j, i, '');%[num2str(i) ',' num2str(j)]);
                    set(stringsh(i,j), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
            end        
            indsA = reshape(1:nr*nR, matDims);
        end
    end
    
    frmIdx_r = nan(1,nFramesPerExtFrame);  % use (fast) indexing to convert frame #s to matrix indices.
    frmIdx_R = nan(1,nFramesPerExtFrame);  % (use nans to ensure that out of bounds entries are not called. [if called, will generate error.])
    frmIdx_r(r_framesToConsider) = 1:length(r_framesToConsider);
    frmIdx_R(R_framesToConsider) = 1:length(R_framesToConsider);
    
    mtx_r_Ind = @(stimId, frmId)   ((stimId-1) * nrFramesConsidered) + frmIdx_r(frmId);
    mtx_R_Ind = @(stimId, frmId)   ((stimId-1) * nRFramesConsidered) + frmIdx_R(frmId);
    
    tic;
    for stimI = 1:nStims
        frameIndsForStimI = find(frameStimIds == stimI);
        
        for frameJ_idx = 1:length(R_framesToConsider)
            frameJ = R_framesToConsider(frameJ_idx);            
            relevantFrameDistsForFrameJ = fliplr(frameJ - r_framesToConsider);
            
            for frameDist_k = 1:length(relevantFrameDistsForFrameJ)
                nFramesAway = relevantFrameDistsForFrameJ(frameDist_k);
                relevantFrameOfKthStim = -nFramesAway + frameJ;
                
                frameIdsForStimK = (frameIndsForStimI + nFramesAway);
                frameIdsForStimK = frameIdsForStimK( ibetween( frameIdsForStimK, 1, nTotalFrames) );

                allStimIdsK = frameStimIds(frameIdsForStimK);
                allStimIdsK = allStimIdsK(allStimIdsK > 0); % some elements can be 0 - if not considering all stimuli (set their values to 0)  
                if ~isempty(allStimIdsK)                    
                    [stimIdsK, stimCount] = uniqueCount(allStimIdsK);

                    rowInds = mtx_R_Ind(stimI, frameJ);
                    colInd  = mtx_r_Ind(stimIdsK, relevantFrameOfKthStim); 
                    A( rowInds , colInd ) = stimCount;
                
                end
            end

        end

    end
    toc;
    
    if showWorking
        set(hM, 'Cdata', A);
        if showDetails
            arrayfun( @(i) set(stringsh(i), 'string', iff(A(i), num2str(A(i)), '') ), indsA);
        end
        3;
    end
    
    if isinteger(A)
        tic;
        [r,c,v] = find(A);
        A = sparse(r,c,double(v), nR, nr);
        toc;
    end
    
    if (nr == nR)
        assert( nnz (A-A') == 0 ); % ie. is symmetric
    end
    
    A = A / length(frameIndsForStimI);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %    decorrelatePSTHs_helper(R, nBinsPerFrame, r_framesToConsider, A)
    [nBinsPerExtFrame, nStims] = size(R);
%     nFramesPerExtFrame = nBinsPerExtFrame/nBinsPerFrame;
    
    R_framesToConsider = 1:nFramesPerExtFrame;
    R_binsConsidered = crossOp([1:nBinsPerFrame]', @plus, [R_framesToConsider-1]*nBinsPerFrame );    
    r_binsConsidered = crossOp([1:nBinsPerFrame]', @plus, [r_framesToConsider-1]*nBinsPerFrame );
    
    R_shapedForMtx = reshape(R(R_binsConsidered(:),:), [nBinsPerFrame, length(R_framesToConsider)*nStims])';

    % Do the actual decorrlation:
    disp('mldivide: (beginning matrix inversion ... wish me luck!)');
    pause(1);
    tic;
    r_shapedForMtx = A \ R_shapedForMtx;
    toc;    
    
    disp('linsolve: (beginning linear solving ... wish me luck!)');
    pause(1);
    tic;
    r_shapedForMtx2 = linsolve(A, R_shapedForMtx);
    toc;
    
    
    r  = reshape(r_shapedForMtx',  [nBinsPerFrame*length(r_framesToConsider), nStims]);
    r2 = reshape(r_shapedForMtx2', [nBinsPerFrame*length(r_framesToConsider), nStims]);
    save solvedMatrices_inv.mat r r2;
    
    % optionally: put r back into matrix the same size as the original R:
    r_big = zeros(size(R));
    r_big(r_binsConsidered(:),:) = r;
    r = r_big;



end


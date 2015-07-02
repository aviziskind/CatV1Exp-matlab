function [nFramesPerPres_Runs, endStatus, runs, presOK] = dbParseSyncs(Did, syncs, dbugFlag)
    % Parse the syncs and try to determine which presentations were shown, 
    % and how many frames were in them.

    % endStatus: 0 if matched syncs to frames perfectly
    %            1 if matched syncs to frames after some fixing
    %            2 if was able to match number of presentations, but not # of frames within each presentation.
    %            -1 if unable to match number of presentations.
        
    
%     Gid = dbLookup('Gid',  'Did', Did);
%     plotSyncsAndSpikes(Gid(1));
    hnd = dbOpenExpDb;
    dbug = exist('dbugFlag', 'var') && ~isempty(dbugFlag);

    runs = [];
    nFramesPerPres_Runs = [];
    endStatus = -1;
    
    % 1. Get the Syncs, and DB fields
    tableName = getDatabaseTableForDid(Did);
    if ~exist('syncs', 'var') || isempty(syncs);    
        syncs = dbGetSyncs('Did', Did, 'tick', 1, 1);
    end
    syncs = syncs(:);
    if dbug
        [stimulusType subType] = getStimulusTypeForDid(Did);
        fprintf(['[' stimulusType ':' subType ']']);
    else
        stimulusType = getStimulusTypeForDid(Did);
    end
        
    if isempty(syncs)
        fprintf('bad syncs:  exiting...');
        return;
    end
    
    nSustFramesPerPres_DB = getFieldsFromDatabaseTable(dbOpenExpDb, 'LNG_N_SUSTAINED_FRM', tableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
    if dbDoesFieldExist(hnd, 'LNG_N_FADEIN_FRM', tableName)
        [nFadeIn, nFadeOut] = getFieldsFromDatabaseTable(dbOpenExpDb, {'LNG_N_FADEIN_FRM', 'LNG_N_FADEOUT_FRM'}, tableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
        nSustFramesPerPres_DB = nSustFramesPerPres_DB + nFadeIn + nFadeOut;
    end
    
    numOrigDBpres = length(nSustFramesPerPres_DB);
    presOK = true(1,numOrigDBpres);
    unknownNumPres = false;
    
    if ~(any(syncs(1:5) < 10))   % add beginning sync.
        syncs = [0; syncs]; 
    end
    frameLength_ticks = getFrameLength('Did', Did, 'tick');

    % Create first attempt at 'runs' variable
    dSyncs_frm = round ( diff(syncs) / frameLength_ticks );
    dSyncs_frm(dSyncs_frm == 0) = 1;  %if any really short frames - round up to1

    isFrame = double( dSyncs_frm <= 1);
    isFrame(~isFrame) = dSyncs_frm(~isFrame);  % if > 1 frame, put in exactly how many framelengths it is.
    runs = runLengths(isFrame);

    runs0 = runs; %#ok<NASGU>

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% Helper functions: %%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%% MAIN PART: FIXING 'RUNS' & JOINING & SEPARATING  PRESENTATIONS

    function ns = get_nFramesPerPres_Runs
         indPres = runs(:,1) == 1; 
         ns = runs(indPres,2);
    end

    function removeRows(rs)        
        runs = runs( setdiff(1:size(runs,1), rs), : ) ; 
    end

    function joinWithNextPres(i)
        assert( mod(i, 1) == 0 ); % whole number
        % 1. Find out how many frames are in the (i+1)th presentation.
        indGapAndNextPres = [(i)*2+1, (i+1)*2]; 
        nFramesInGapAndNextPres = sum( runs(indGapAndNextPres, 2) );        

        % 2. Add those frames to the current one (together with inter-presentation frames).
        runs( i*2, 2 ) = runs( i*2, 2 ) + nFramesInGapAndNextPres;

        % 3. Delete the next presentation.        
        removeRows( i*2 + [1,2] );        
    end
    

    allAreAlmostEqual = @(x,y, n) (length(x) == length(y)) && all(abs(x - y) <= n);
    function tf = almostEqlExceptForAFew(frms1, frms2)
        tf = false;
        if length(frms1) ~= length(frms2)            
            return;
        end
        isMatch = (frms1 == frms2);
        discrepSign = sign(frms1-frms2);
        anyMisassignments = any( abs(diff(discrepSign)) == 2);  % ie. any 'too short' followed by 'too long', or visa versa.
        if (nnz(isMatch)/length(isMatch) > .95) && ...
                ~anyMisassignments  % ie. all discrepancies are local, and the mistakes are most likely not due to parsing errors.
            tf = true;
        end        
    end

    function tf = matchesWell(nFramesPerPres)
        tf = isequal(nFramesPerPres, nSustFramesPerPres_DB ) || ...  % Ideal situation
               all([length(nSustFramesPerPres_DB) length(nFramesPerPres)] == 1)  || ...      % these are always easy to match up.
                ( strcmp(stimulusType, 'Grating') && ...         % gratings are continuous showings of gratings with the same params, so ok if off by a few frames here and there.
                     (   allAreAlmostEqual(nFramesPerPres_Runs, nSustFramesPerPres_DB, 3) ) ... 
                      || almostEqlExceptForAFew(nFramesPerPres_Runs, nSustFramesPerPres_DB) );                      
    end

    function tf = matchesFrameNum(nFramesPerPres)
        tf = (length(nFramesPerPres) == length(nSustFramesPerPres_DB));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%% FIXING 'RUNS' & JOINING & SEPARATING  PRESENTATIONS %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Step 0. Get meanNumSustFrames
    nFramesPerPres_Runs = get_nFramesPerPres_Runs;    
    if ~(all(nSustFramesPerPres_DB) == 0)
        meanNumSustFrames = median(nSustFramesPerPres_DB);
    else
        meanNumSustFrames = median(nFramesPerPres_Runs);
        if mod(meanNumSustFrames, 10)
            if length(nSustFramesPerPres_DB) == 1                
                meanNumSustFrames = sum(nFramesPerPres_Runs);
            else
                keyboard;
            end
        end
    end
    

    % STEP 1. Remove double gaps    
    while true
%         interPresGaps = (runs(:,1) > 1) & (runs(:,2) == 1);  %ie. "one occurence of extended frame".
        interPresGaps = (runs(:,1) > 1);  %ie. "N occurences of extended frames".
        doubleGaps = interPresGaps(1:end-1) & interPresGaps(2:end); % two consecutive 'gaps' (points to first occurence).
        doubleGapInds = find(doubleGaps);  % indices of first gap followed by a second one.
        doubleGapVals = runs(doubleGapInds, 1);  % the length (in time) of the gap.
        %remove double gaps:
        if isempty(doubleGapInds)
            break;
        else  %eg Did = 191            
%             gapPos = sort(doubleGapInds-1, 'descend');  
            gapPos_idx = indmin(doubleGapVals);  
            gapPos     = doubleGapInds(gapPos_idx);
            % more important factor: we don't want to add frames to a presentation that makes it > meanNumSustFrames
            % other factor: we want to join the smallest fragments together
            for g = gapPos(:)';
                atEndOfArray_Factor = [0 0];
                aboveMeans_Factor = [0 0];
                smallerFrags_Factor = [0 0];                
                atEndOfArray_Factor(1) = 4*(g == 1);
                atEndOfArray_Factor(2) = 4*(g >= size(runs,1)-1);
                aboveMeans_Factor(1) = 2*((g == 1)              || (sum( runs([g,   g-1], 2) ) > meanNumSustFrames));
                aboveMeans_Factor(2) = 2*((g >= size(runs,1)-1) || (sum( runs([g+1, g+2], 2) ) > meanNumSustFrames));                
                smallerGapInd = indmin( runs([g, g+1],1) );
                smallerFrags_Factor(smallerGapInd) = smallerFrags_Factor(smallerGapInd)+1;
                Q = atEndOfArray_Factor + aboveMeans_Factor + smallerFrags_Factor;
                
                if Q(1) < Q(2) 
                    runs(g-1,2) = runs(g-1,2) + runs(g, 2);
                    removeRows(g);                
                else
                    runs(g+2,2) = runs(g+2,2) + runs(g, 2);
                    removeRows(g);
                end            
            end
        end
    end
    
%     runs_noDoubleGaps = runs;    
            
    nPres_syncs = length(nFramesPerPres_Runs);
    nPres_DB = length(nSustFramesPerPres_DB);
    
    % STEP 3a:  Check : If nPres from syncs < nPres from DB : likely that experiment was was interrupted  
    if (nPres_syncs < nPres_DB) && ...
            (sum(nFramesPerPres_Runs) < sum(nSustFramesPerPres_DB))
        if dbug
            fprintf('nPres_syncs (%d) < nPres_DB (%d) (exp interrupted?)', length(nFramesPerPres_Runs), length(nSustFramesPerPres_DB));
        end
        unknownNumPres = true;        
        numShownPres = length(nFramesPerPres_Runs);
        nSustFramesPerPres_DB(numShownPres+1:end) = [];
        % known bad ones:
        if any(Did == [2913])
            return;
        end
        
        if dbug
            disp(runs);
        end
    end
    
    %  STEP 3b:  Check : If nPres from DB < nPres from syncs
    %              (AND  nFrames_DB < nFrames_Syncs) :  likely that some presentations were left out of the database.

    if (nPres_DB < nPres_syncs) && ...
        ( (sum(nSustFramesPerPres_DB) < sum(nFramesPerPres_Runs) -length(nSustFramesPerPres_DB) - 10)  );
        % [[[ eg Did = 1622 (only 1/2), 2639 (only(3/4)), [* There is a script for adding these to the database* ] ]]]  
        % [[[ eg Did = 474, 479, 480 [* Ticks repeated. Fix added to dbGetSyncs to handle these cases .*] ]]]
        % eg Did = 76:  just have extra frames in 1 pres. but otherwise ok.
        abortTry = true;        
        if any(Did == [76])  
            joinWithNextPres(156);
            joinWithNextPres(156);
            nFramesPerPres_Runs = get_nFramesPerPres_Runs;
            abortTry = false;
        end        
        if abortTry            
            if dbug
                fprintf('It seems some presentations are missing from the database.'); 
            end
%             keyboard;
            endStatus = -1;
            return;        
        end
    end     
        

    %%% STEP 4: JOIN FRAGMENTED PRESENTATIONS
    
    if matchesWell(nFramesPerPres_Runs)
        if dbug
            disp(runs);
            fprintf('[Matches database]');
        end
%         fprintf('[Matches database !]');
        endStatus = 0;    
    else
            
        if dbug
            fprintf('presentations broken up. Attempting to join them ...\n');
            disp('Before:');
            disp(runs);
        end        
        % The following algorithm joins fragmented presentations by
        % (1) finding the  presentations which have the smallest (false) 
        % inter-presentation gaps between them, but whose total sum is less 
        % than the meanNumSustFrames. If none can be found, then we:
        % (2) join consecutive presentations whose sum is as small as 
        % possible above meanNumSustFrames.
                       
        preBlank_tmp = runs(1,1);
        runs(1,1) = 10*max(runs(:));  % hack to make sure first gap isn't removed.
              
        
        % step 1. join presentations, as long as they don't add up to more
        % than the expected number of frames per presentation 
        nFramesPerPres_Runs = get_nFramesPerPres_Runs;        
        sumOfConsecutivePres = [nFramesPerPres_Runs(1:end-1) + nFramesPerPres_Runs(2:end)];                
        
        nPresOrig = length(nFramesPerPres_Runs);
        joinedPresList = num2cell(1:nPresOrig);
        
        while ((length(nFramesPerPres_Runs) > length(nSustFramesPerPres_DB)) || unknownNumPres) ...
                || any(sumOfConsecutivePres < meanNumSustFrames);
            
            gapInds =  runs(:,1) > 1;
            gaps = runs(gapInds, 1);
            [tmp, smallestInds] = sort(gaps);
            i = 1;
            foundPresToJoin = false;
            
            % try (1): finding smallest gaps.
            while ~foundPresToJoin && (i <= length(smallestInds))
                indPres = smallestInds(i);
                if (indPres < length(nFramesPerPres_Runs)) ...
                    && ((nFramesPerPres_Runs(indPres) + nFramesPerPres_Runs(indPres+1) <= meanNumSustFrames))
                    foundPresToJoin = true;
                end
                i = i+1;
            end
            
            % (2): finding consecutive presentations with smallest excess.
            if ~foundPresToJoin
                indPres = indmin(sumOfConsecutivePres);    
            end
            
            joinWithNextPres(indPres);
            joinedPresList{indPres} = [joinedPresList{[indPres, indPres+1]}];
            joinedPresList(indPres+1) = [];
            assert(isequal([joinedPresList{:}], 1:nPresOrig));
            
            nFramesPerPres_Runs = get_nFramesPerPres_Runs;
            sumOfConsecutivePres = [nFramesPerPres_Runs(1:end-1) + nFramesPerPres_Runs(2:end)];
                
        end
        
        runs(1,1) = preBlank_tmp;
        dispRuns = runs;
        presOK = cellfun(@(x) length(x == 1), joinedPresList);
        
        if dbug
            disp('After:');
            disp(dispRuns);
        end
                
        
        if matchesWell(nFramesPerPres_Runs)
            endStatus = 1; % ok after fix #1        
    %         fprintf('[Method 1 worked !]');        
    
        elseif matchesFrameNum(nFramesPerPres_Runs)  % couldn't match presentations exactly, but at least matched # of presenations
            endStatus = 2;
        
        
        end
    end
    
    
    if length(nFramesPerPres_Runs) < numOrigDBpres  % b/c experiment was interrupted
        nFramesPerPres_Runs(end+1:numOrigDBpres) = 0; % add back original presentations (with 0 frames).
        if endStatus == 0
            endStatus = -1;
        end
    end
    
    
end





% function joinedPres = findBestPresJoining(splitPres, nPres, nEachPres)
%     splitPres = [100   50 40 10   100  50 50  40 60  100 30 40 28 30];
%     nPres = 9;
%     nEachPres = 100;
% 
% 
%     joinedPres = splitPres;
%     indJoinedPres = cell(1,nPres);
%     ind = 1;
%     
%     fullPres = find(splitPres == nEachPres);
%     splitGroups = continuousGroupings(1:length(splitPres), fullPres);
%     
%     while length(joinedPres) < nPres
%         
%         
%         
%     end
%     
% 
% 
% 
% 
% 
% 
%     
% 
% end
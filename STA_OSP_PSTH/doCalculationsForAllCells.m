function doCalculationsForAllCells(groupIdsToDo, cellIdsToDo)
    global psthStatsSettings
    if isempty(psthStatsSettings)
        setPsthGlobalSettings;
    end

	acceptableErrorIds = {'db:noSpikes', 'db:noSid', 'db:noSyncs', 'db:badSyncs', 'db:invalidTbTe', ...
                          'db:nFramesTxtFileMismatch', 'db:syncsTbTeMismatch', 'db:syncNframesMismatch', ...
                          'stim:interrupted', 'stim:singleTrial', 'db:fps', 'db:SinglePhaseMovie', 'db:corruptedSpikeData', 'db:flawedSyncs'}; % 

    doProgressBar = 0;
                      
    if nargin < 1
%     groupIdsToDo = [337 264 282 1048 1141 338 265 283 481 463 432 465 466 473 474 476 483 484 488 490 492 493 494 495 496 497 502 503 504 505 506 507 511 514 515 516 517 520 521 522 523 525 526 528 530 529 531 537 538 543 544 545 546 453 454 549 550 551 552 553 554 558 559 560 563 565 566 567 568 569 573 574 575 576 598 599];
%     groupIdsToDo = [849 852 853 855 856 873 917 918 920 924 928 929 931 932 933 1008 1009 1012 1013 1105 1112 1129 1141 1143 1150 1156 1163 1164 1438 1512 1519];
%     groupIdsToDo = [1048];
%     groupIdsToDo = [461];
%     groupIdsToDo = [4470];
%     groupIdsToDo = [2460 3073 4344 4362 4768 5100 5146]; % groups with sites not 100% ok.
%     groupIdsToDo = [1857 1875 1877 2268 2270 2272 2288 2290 2292 2302 2304 2312 2314 2344 2346 2352 2354 2362 2368 2374 2376 2416 2422 2432 2434 2440 2448 2458 2462 2464 2478 2480 2484 2486 2500 2502 2504 2506 2514 2516 2518 2520 2542 2544 2546 2548 4460 4462 4466 4468 4470 4474 4476 4480 4482 4486 4488 4492 4494 4498 4502 4506 4508 4514  4518 4520 4522 4538 4544 4546 4548 4590 4596 4702 4706 4708 4712 4714 4716 4718 4722 4726 4730 4732 4738 4744 4746 4778 4784 4788 4794 4796 4798 4878 4882 4888 4890 4904 4910 4930 4936 5084 5086 5112 5114];    

%     groupIdsToDo = [1877 2312 2346 2761 2771 4494 4538];
%     groupIdsToDo = [2771 4494 4470];
%         groupIdsToDo = [1648];
    end

    global nCellsCompleted nFramesCmp;
    nFramesCmp = [];   
    if exist('stop_program', 'file')
        movefile('stop_program', 'stop_progra')
    end

    clustGrp = curGroupingType('');  % 'cells' or 'clusters'    
    
    stimTypes = {'movie', 'noise', 'grating', 'mseq'};
    typesToDo = [1];

    % for specific stimulus types:
    restrictToFlashGratingMovies = true;
    restrictToDrifting_SpfBatch = false;
    restrictToDrifting_OriBatch = false;    
    restrictToDriftingGratingWithMultipleSpf = false;
    restrictToDriftingGratingWithSingleTpf = true;    
    restrictToDriftingGratingWithNPhases = [];%40, 60];
    ignoreSingleGratingStims = true;
    
    
    
    if (nargin > 0) 
        allSD = siteDataFor(groupIdsToDo);
        L = 20;
        grpStimTypes = arrayfun(@(s) s.stimType(1:L), allSD, 'un', 0);
        uType = unique( grpStimTypes );
        if length(uType) > 1
            error('Multiple types');
        end
        if strncmp(uType, 'Movie:Flashed_Gratings', L);
            typesToDo = 1;
        elseif strncmp(uType, 'Grating:Spatial Frequency Batch', L);
            typesToDo = 3;
            restrictToDrifting_SpfBatch = 1;
        elseif strncmp(uType, 'Grating:Orientation Batch', L);
            typesToDo = 3;
            restrictToDrifting_OriBatch = 1;
        end
    end
    

    % instructions:
    justValidateSites = false;
    
    % how to choose which sites:
    whichCellsToDo = 'GidList'; % options: 'all', 'notdone', 'failed', 'notfailed', 'GidList';
    useOldDataIfExists = true;
    onlyDoOneCellFromEachGroup = false;        
    
%     acceptUnknownErrors = false;
      
    nCellsCompleted = 0;
    runsBetweenSaves = 500;
    
    
    cellFunction = iff(justValidateSites, @(G, c, o) testSiteValidity(G, 'err'), @calculatePSTH_STAs_OSP_ForOneCell);    
    varnamePrefix = iff(justValidateSites, 'siteTest', 'celldata');
    if justValidateSites
        onlyDoOneCellFromEachGroup = true;        
        isFailedCell = @(celldata) isstruct(celldata);
    else
        isFailedCell = @(celldata) isfield(celldata, 'id') || (isfield(celldata, 'OSP') && ~isfield(celldata.OSP, 'R_full'));        
    end
    
    groupIdsProvided = (nargin >= 1) && ~isempty(groupIdsToDo);    
    cellIdsProvided = (nargin >= 2);
    if groupIdsProvided && (length(groupIdsToDo) == 1) && cellIdsProvided
        groupIdsToDo = groupIdsToDo(ones(1, length(cellIdsToDo)));
    end
    if cellIdsProvided && (length(groupIdsToDo) ~= length(cellIdsToDo))
        error('cellIds provided must match length of groupIds (unless groupIds is of length 1)')
    end    
            
    for ti = typesToDo
        stimType = stimTypes{ti};
        if length(typesToDo) > 1
            disp([' ---------------------------------------------- ']);
            disp([' NOW DOING ' upper(stimType) ' ' upper(clustGrp) ' ... ']);
            disp([' ---------------------------------------------- ']);
        end
        
        ext = '';        
        if strcmp(stimType, 'movie') && restrictToFlashGratingMovies
            ext = '_fg';
        end
        if strcmp(stimType, 'grating') && restrictToDrifting_SpfBatch
            ext = '_dSf';
        elseif strcmp(stimType, 'grating') && restrictToDrifting_OriBatch
            ext = '_dOr';
        end
        
        
        basepath = CatV1Path;
        if justValidateSites
            indivCellRecords_filename = ['testResults_' stimType clustGrp ext '_Sites.mat'];
        else
            [basepath, indivCellRecords_filename] = getFileName('indiv', [stimType ext]);            
%             clustGrp1 = clustGrp; clustGrp1(1) = upper(clustGrp1(1));
%             indivCellRecords_filename = ['indiv' clustGrp1 curCellsType '_' stimType ext '.mat'];
        end                
        
        if exist([basepath indivCellRecords_filename], 'file')
            fprintf('[[ LOADING ... '); tic; load([basepath indivCellRecords_filename]); fprintf('done ]] '); toc; 
        end

        % retrieve master list of Gids & cellIds if explicit pairs not provided
        if ~cellIdsProvided
        
            cellGroupsFile = getFileName('Groups', [stimType ext]);
            S = load(cellGroupsFile);
            groupData = S.([stimType 'Groups' ext]);        
    %         eval([stimType 'Groups' ext ' = groupData;']);

            if strcmp(stimType, 'grating')  && (nargin < 1)
                if ignoreSingleGratingStims
                    inds_notSingle = findInStructArray(groupData, 'stimType', [], @(s) ~strncmp(s, 'Grating:Single', 14));
                    groupData = groupData(inds_notSingle);
                end        
                if restrictToDriftingGratingWithMultipleSpf 
                    ind_multSpf = findInStructArray(groupData, 'spPeriod_pix', [], @(x) length(x) > 1);    
                    groupData = groupData(ind_multSpf);
                end        
                if restrictToDriftingGratingWithSingleTpf            
                    ind_singleTpf = findInStructArray(groupData, 'tempPeriod_sec', [], @(x) length(x) == 1);
                    groupData = groupData(ind_singleTpf);
                end            
                if ~isempty(restrictToDriftingGratingWithNPhases)
                    3;
                    ind_GoodNphases = findInStructArray(groupData, 'spPh_deg', [], @(x) any(length(x) == restrictToDriftingGratingWithNPhases));
                    groupData = groupData(ind_GoodNphases);
                end

            end
        end
        
        if (length(groupIdsToDo) == 1) && (length(groupIdsToDo) ~= length(cellIdsToDo))
            groupIdsToDo = groupIdsToDo * ones(1, length(cellIdsToDo));            
        end
        
        if ~cellIdsProvided || (length(groupIdsToDo) ~= length(cellIdsToDo))                    
            if strcmp(whichCellsToDo, 'GidList')
                idx = arrayfun(@(Gid) any(groupIdsToDo == Gid), [groupData.Gid]);
                groupData = groupData(idx);
            end

            allStimGids = num2cell([groupData(:).Gid]');
            if cellIdsProvided
                allStimCellIds = {double(cellIdsToDo)};
            else            
                if ~onlyDoOneCellFromEachGroup
                    allStimCellIds = arrayfun(@(grp) [double(grp.cellIds) 100], groupData, 'un', 0);
                else
                    allStimCellIds = arrayfun(@(grp) double(grp.cellIds(1)), groupData, 'un', 0);
                end
            end        
            allCells = cellfun(@(grp, cells) [grp*ones(1, length(cells)); cells], allStimGids, allStimCellIds, 'un', 0);
            allCells = [allCells{:}]';

        else
            allCells = [groupIdsToDo(:), cellIdsToDo(:)];
        end
        
        s = who([varnamePrefix '_*']);
        [GidsDone cellIdsDone] = cellfun(@getGidCellIdFromVarname, s);
        allCellsDone = [GidsDone(:), cellIdsDone(:)];
        
        if any(strcmp(whichCellsToDo, {'failed', 'notfailed'}))
            s_details = whos([varnamePrefix '_*']);
            failedCellInds = false(1,length(s));
            for i = 1:length(s)
                failedCellInds(i) = isFailedCell(eval(s{i}));
                if failedCellInds(i)
                    3;
                end
            end
            [failedGids, failedCellIds] = cellfun(@getGidCellIdFromVarname, s(failedCellInds));
            allFailedCells = [failedGids(:), failedCellIds(:)];
        end
        
        switch whichCellsToDo
        
            case {'all', 'GidList'}
                allCellsToDo = allCells;    
                allCellsDone = [];
                
            case 'notdone'
                allCellsToDo = setdiff(allCells, allCellsDone, 'rows');

            case 'failed',
                allCellsToDo = allFailedCells;

            case 'notfailed',
                allCellsToDo = setdiff(allCells, allFailedCells, 'rows');
                                
        end
           
        onlyDoMU = 0;
        if onlyDoMU
           idx_mu = (allCellsToDo(:,2) <= 0) | (allCellsToDo(:,2) == 100);
           allCellsToDo = allCellsToDo(idx_mu,:);
            
        end
        
        % calculate how many cells are already done, & how many we need to do        
        nCellsDone = size(allCellsDone,1);
        nCellsToDo = size(allCellsToDo,1);
        
        CellStartAt = 1;
        if doProgressBar
            progressBar('init=', [nCellsDone, nCellsDone+nCellsToDo]);
        end
        for cell_i = CellStartAt:nCellsToDo
            if doProgressBar
                progressBar;
            end
            [Gid, cellId] = deal(allCellsToDo(cell_i,1), allCellsToDo(cell_i,2));

%             group_idx = find([groupData.Gid] == Gid, 1);
%             if isempty(group_idx)
%                 continue;
%             end
%             thisGroupData = groupData(group_idx);

            
            varname = getName(varnamePrefix, Gid, cellId);

            oldCellData = [];
            if exist(varname, 'var') && useOldDataIfExists
                oldCellData = eval(varname);
            end
            tic;
            if   any(strcmp(whichCellsToDo, {'all', 'GidList'})) || ...
                (strcmp(whichCellsToDo, 'notdone') && isempty(oldCellData)) || ...                
                (strcmp(whichCellsToDo, 'failed') && (isempty(oldCellData) || isFailedCell(oldCellData))) || ...
                (strcmp(whichCellsToDo, 'notfailed') && ~isFailedCell(oldCellData));
%                 try
                    oktxt = '_';%iff(all(thisGroupData.presOK), '-', '!');
                    cellId_str = iff(cellId == 100, '*', num2str(cellId));
                    fprintf(['[' outOf(cell_i, nCellsToDo) '] {%d:%s}  [' oktxt ']'], Gid, cellId_str);
                    newCellData = cellFunction(Gid, cellId, oldCellData);

                    eval([varname ' = newCellData;']);
                    cell_accepted = iff(newCellData.OSP.stats.isRep, '=', 'x');
                    disp([ '[' cell_accepted  '] ' num2str(toc, '%.2f') ' seconds.']);
%                  catch Merr                        
%                      fprintf(' ... failed!');                        
%                      if any( strcmp(Merr.identifier, acceptableErrorIds) ) || acceptUnknownErrors % i.e. is a known error.
%                          newCellData = struct('id', Merr.identifier, 'txt', Merr.message);
%                          disp(Merr.message);
%                          eval([varname ' = newCellData;']);
% %                      else
% %                          rethrow(Merr);
%                      end
%                 end
            else                
                disp(['[' outOf(cell_i, nCellsToDo) '] : Skipping this cell:   ']);
            end

            if (mod(nCellsCompleted+1, runsBetweenSaves) == 0)
                %%
                fprintf('[[ SAVING ... '); tic; save([basepath indivCellRecords_filename], [varnamePrefix '_*'], '-v6'); fprintf('done ]]'); toc; 
                nCellsCompleted = nCellsCompleted + 1;
            end

            if exist('stop_program', 'file');
%                 cleanupAndFixVariables;
                fprintf('[[ SAVING ... '); tic; save([basepath indivCellRecords_filename], [varnamePrefix '_*'], '-v6'); fprintf('done ]]'); toc; 
                disp(['Completed ' num2str(nCellsCompleted)  ' ' clustGrp ])
                return;
            end
                
        end  % of loop for cell                  
        if doProgressBar
            progressBar('done');
        end

%         cleanupAndFixVariables;
        fprintf('[[ SAVING results to file %s ... ', indivCellRecords_filename); tic; save([basepath indivCellRecords_filename], [varnamePrefix '_*'], '-v6'); fprintf('done ]]'); toc;
        
    end % of loop for each stimulus type
                       
    disp(['Completed ' num2str(nCellsCompleted) ' ' clustGrp ])

    dbGetCellSpkStimHists('save');
    getPSTHwindowData('save');
    getOspDataForPsthWindow('save');
    
    if justValidateSites
        A = nFramesCmp
    end
    
    if 0
       save([basepath indivCellRecords_filename], [varnamePrefix '_*'], '-v6'); 
        dbGetCellSpkStimHists('save');
        getPSTHwindowData('save');
        getOspDataForPsthWindow('save');
        
    end
    
end


% Gids with fgcomp movie: [ 2625 2665 2677 2723 2737 2747 2775 2791 2805 2887 2921 2957 2981 3047 3061 3077 3091]



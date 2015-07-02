function S_amend = dbGetSyncAmendments(idType, idVal)
%     S_amend = dbGetSyncAmendments(idType, idVal)        
%     S_amend = dbGetSyncAmendments(Gid)

    persistent allSyncAmendments saveCount 
    
    syncAmendments_file = [CatV1Path 'MatLabDB_avi' filesep 'allSyncAmendments.mat'];
        
    redo_all = 0;
    redo_current = 0;
    saveCountSpacing = 50;

    if strcmp(idType, 'save')
        save(syncAmendments_file, 'allSyncAmendments', '-v6');                
        saveCount = 0;
        return;
    end        

    if nargin == 1 ; %
        idVal = idType; 
        idType = 'Gid';
    end    

    
    if isempty(allSyncAmendments)        
        if exist(syncAmendments_file, 'file') && ~redo_all
            S_file = load(syncAmendments_file);
            allSyncAmendments = S_file.allSyncAmendments;
        else
            allSyncAmendments = struct;
        end        
        saveCount = 0;        
    end
        
    Gid = dbLookup('Gid', idType, idVal);    
    assert(length(Gid) == 1);
    if flashedOrDrifting(Gid)
        sync_fld_name = sprintf('SyncAmendments_Gid_%d', Gid);        

        if (~isfield(allSyncAmendments, sync_fld_name) || redo_current)
            S_amend = calcSyncAmendments(Gid);
            if ~isempty(fieldnames(S_amend))
                3;
            end

            allSyncAmendments.(sync_fld_name) = S_amend;
            saveCount = saveCount + 1;

            if saveCount > saveCountSpacing
                save(syncAmendments_file, 'allSyncAmendments', '-v6');        
                saveCount = 0;
            end                
        end

        S_amend = allSyncAmendments.(sync_fld_name);
    else
        S_amend = struct;
    end
    

end

function S_amend = calcSyncAmendments(Gid)
    S_amend = struct;

    sd = siteDataFor(Gid);
    Did = sd.Did;
    nPres = length(sd.presIds);
    frameLength_ms = sd.frameLength_ms;
    frameLength_tick = dbConvertTimeMeasures(Did, frameLength_ms, 'ms', 'tick');

    retrieveWithoutAmendments_flag = 1;
    syncTimes_ticks = dbGetSyncs('Gid', Gid, 'tick', [], retrieveWithoutAmendments_flag);

    if ~issorted(syncTimes_ticks)
        % for Gid = [437, 438, 439], (Sid = [376, 381, 382]),  for some reason, the syncs are
        % repeated twice  (so just take the first half of the syncs)
        assert(any(Gid == [437, 438, 439]));
        S_amend.cutInHalf = true;
        
        return;
    end
        
    
    si = getGratingStimType(Gid);
    nFramesPerPres_correct = (si.nOri * si.nSpf * si.nPh)*si.nTrials / nPres;
        
    [presStartSyncIds, presEndSyncIds] = getTbTeIdxFromTicks(syncTimes_ticks, frameLength_tick);
    
    if length(presStartSyncIds) ~= nPres
        error('Incorrect # of presentations');
    end
        
    nFramesEachPres_displayed = presEndSyncIds-presStartSyncIds;

    
    nFrames_discrepancy = (nFramesEachPres_displayed(:)' - nFramesPerPres_correct);
    idx_discrepancy = find(nFrames_discrepancy);                       
        
    if ~isempty(idx_discrepancy)
        syncsToRemove = cell(1, nPres);
        syncsToAdd    = cell(1, nPres);

        for pres_idx = idx_discrepancy
            nDiscrep = nFrames_discrepancy(pres_idx);
            if nDiscrep > 0  % too many frames: remove the extraneous frames.                    
                if (Gid == 513) 
                    % 2 electronic blips in experiment which caused spurious syncs.
%                         idx = find(dt < frameLength_tick/2)
%                         presEndSyncIds(pres_idx) - [0:discrep-1];                        
                    syncsToRemove{pres_idx} = [6377004 6572835];

                else
                    error('Not yet implemented');
                end

            elseif nDiscrep < 0 % too few frames: fill in the missing frames:

                pres_lastTick = syncTimes_ticks (presEndSyncIds (pres_idx) );
                syncsToAdd{pres_idx} = pres_lastTick + frameLength_tick * [1:abs(nDiscrep)];

            end

        end
        
        syncsToRemove = unique([syncsToRemove{:}]);        
        if ~isempty(syncsToRemove)
            S_amend.syncsToRemove =syncsToRemove;
        end
        
        syncsToAdd = [syncsToAdd{:}];
        if ~isempty(syncsToAdd)
            S_amend.syncsToAdd =syncsToAdd;
        end

        show = 0;
        if show
            figure(Gid); clf;
            plot(syncTimes_ticks, ones(size(syncTimes_ticks)), '.') 
            drawVerticalLine(syncTimes_ticks(presStartSyncIds), 'color', 'g')
            drawVerticalLine(syncTimes_ticks(presEndSyncIds), 'color', 'r')

            drawVerticalLine(syncsToAdd, 'color', 'k')
        end
        
        chk = 1;
        if chk
            syncTimes_ticks_fixed = sort([syncTimes_ticks; syncsToAdd(:)]);
            syncTimes_ticks_fixed = setdiff(syncTimes_ticks_fixed, syncsToRemove);
            [tb_idx, te_idx] = getTbTeIdxFromTicks(syncTimes_ticks_fixed, frameLength_tick);
            nFramesEachPres_fixed = te_idx-tb_idx;
            assert(all( nFramesEachPres_fixed == nFramesPerPres_correct));                        
        end
        
    end
    3;
        
end

function [tb_idx, te_idx] = getTbTeIdxFromTicks(syncTimes_ticks, frameLength_tick)
    dt = diff(syncTimes_ticks(:));
    idx_interval = find(dt > frameLength_tick*2);

    tb_idx = [1; idx_interval+1];
    te_idx = [idx_interval; length(syncTimes_ticks)];
    
end


%         [uDisp, dispCount] = uniqueCount(nFramesEachPres_displayed);        
%         if all(uDisp <= nFramesPerPres_correct)
%             s = 'MISSING';
%         elseif all(uDisp >= nFramesPerPres_correct)
%             s = 'EXTRA';
%         else
%             s = 'MIXED';
%         end
%         fprintf(' Gid = %d: %s : %s \n', Gid, s, sprintf(' %d (n=%d)', [uDisp(:), dispCount(:)]') );
%         return;

%{
allGids = [getAllGids('f'), getAllGids('o'), getAllGids('s')];
progressBar('init-', length(allGids), 60);
for i = 1:length(allGids)
    dbGetSyncAmendments(allGids(i));
    progressBar(i);
end
%}

%{
 errorGids = [513 605 607 608 610 612 615 616 619 620 624 628 631 633 634 640 644 650 652 662 663
 664 665 667 668 670 675 676 678 682 687 698 700 701 703 709 710 713 714 716 719 723 724 728 730 
732 734 735 736 741 742 743 755 757 761 763 765 768 771 772 782 784 788 790  807 808 810 813 814
 817 818 820 826 830 831 833 834 835 837 842 844 851 858 860 861 869 871 875 908 1002 1004 1005 
1101 1103 1107 1108 1109 1110 1111 1114 1115 1116 1118 1119 1121 1122 1124 1125 1126 1130 1131
 1134 1135 1136 1138 1140 1144 1145 1146 1147 1151 1152 1153 1159 1160 1162 1167 1169 1170 1171
 1173 1175 1176 1192 1197 1200 1216 1218 1435 1453 1464 1468 1472 1473 1475 1476 1483 1484 1494
 1498 1500 1511 1513 1516 1518 1520 1527 1528 1530];

%}
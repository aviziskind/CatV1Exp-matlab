function dbFixSustainedDisplayedFrms(lims)

     DidsThatNeedFixing = {[2460 1933 2749 2840 2849 2960 2996 3001 3006 3007 3029 3062 3140 3200 3310], ...
         [ 70 71 111 122 123 131 132 135 136 138 144 155 165 177 183 195 196 200 221 227 230 244 246 251 263 274 287 290 311 312 313 314 315 320 332 333 350 351 372 381 391 454 456 482 818 828 1008 1009 1948] ...
         [ ...191 337 419 467 639 641 649 661 673 677 719 724 729 752 790 792 803 669 694 756  766 812 1007 1034 1051 1072 1075 1087 1100 1121 1092 1011 1126 1128 1202 1212 1220 1404 1411 1427 ...
            ... 1433 1443 1449 1451 1458 1502 1503 1516 1517 1527 1529 1550 1560 1564 1640 1677 1680 1697 1714 1754 1789 1808 1892 1939 1964 2339 2366 2395 2445 2449 2460 2623 2855 2913 2922 2928 3026 3028 3259 3747 3768 ...
            500  420  621  627  629  633  636  644  647  652  656  667  678  683  688  692  705  715  732  735  761  764  767  771  777  781  786  793  796  799  809  820  672  676  680  686   693  708  711  718  721  728  734  738  696  699  702  712  714  740  745  746  751  770  773  780  784  788  798  802  806  823  826  831  834 1016 1017 1019 1023 1024 1029 1031 ...
             1036 1043 1053 1057 1060 1063 1064 1067 1082 1096 1107 1112 1117 1118 1122 1119 1090 1113 1010 1026 1035 1038 1041 1056 1061 1079 1101 1103 1123 1124 1127 1129 1130 1131 1132 1133 1134 1135 1136 1142 1143 1144 1204 1205 1208 1209 1214 1215 1210 1211 1218 1243 1248 1254 1266 1271 1276 1284 1289 1292 1300 1305 1310 1318 1330 1337 1344 1309 1403 1405 1412 1415 ...
            1416 1419 1420 1421 1425 1428 1429 1432 1435 1436 1437 1441 1444 1445 1446 1452 1453 1456 1463 1464 1469 1474 1477 1482 1488 1493 1497 1499 1501 1505 1510 1511 1513 1519 1530 1532 1533 1536 1539 1542 1547 1552 1553 1554 1561 1566 1569 1594 1599 1602 1610 1611 1615 1618 1620 1621 1636 1652 1653 1673 1679 1692 1702 1703 1710 1713 1715 1726 1728 1734 1735 1739 ...
            1743 1747 1748 1750 1751 1758 1759 1761 1765 1767 1769 1773 1775 1783 1786 1787 1788 1791 1793 1794 1795 1798 1802 1803 1805 2445 2460 2913 3026 3028 3274 3659  ], ...
          []};        

      
errorDids = [500 621 627 629 633 636 644 647 652 656 667 672 676 678 680 683 686 688 692 693 696 699 702 705 708 711 712  ...
714 715 718 721 728 732 734 735 738 740 745 746 751 761 764 767 770 771 773 777 780 781 784 786 788 793 796 ...
798 799 802 806 809 820 823 826 831 834 1010 1016 1017 1019 1023 1024 1026 1029 1031 1035 1036 1041 1043 ... 
1053 1057 1060 1063 1064 1067 1082 1090 1096 1103 1112 1113 1204 1205 1208 1403 1405 1415 1416 1419 1420 ...
1421 1428 1429 1432 1435 1436 1437 1441 1444 1445 1446 1453 1456 1463  1464 1469 1474 1477 1493 1497 1499 ...
1501 1510 1511 1513 1530 1532 1533 1547 1552 1553 1554 1561 1566 1569 1594 1599 1602 1618 1620 1728 1710  ...
1739 1743 1759 1750 1747 1748 1802 1803 1805 1751 1758 1769 1773 1775 1786 1788 1791 1793 1795];     
      DidsThatNeedFixing = errorDids;
      
    frameFields = {'LNG_N_DISPLAYED_FRM', 'LNG_N_SUSTAINED_FRM', 'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM', 'LNG_PRESENT_NO'};
    extraGratingFrameFields = {'LNG_N_FADEIN_FRM', 'LNG_N_FADEOUT_FRM'};

    onlyCorrectZeros = false;
    onlyCorrectnDisp = false;
    hnd = dbOpenExpDb;
    
    stimTypes = {'movie', 'noise', 'grating', 'mseq'};
    stimTypesToDo = [3];

    for si = 1:length(stimTypesToDo)
        stim = stimTypesToDo(si);
        stimType = stimTypes{stim};
        stimTable = getDatabaseTableForDid([], stimType);
        disp([' *** STARTING ' upper(stimType) ' STIMULI ***']);
%         Dids = dbGetStimulusDids(stimType);        
    
        Dids = DidsThatNeedFixing; %{stim};
        if nargin < 1
            startAt = 1;
            endAt = length(Dids);
        else
            if length(lims) == 1
                [startAt, endAt] = deal(lims,lims);
            else
                [startAt, endAt] = deal(lims);
            end
        end
        
        progressBar('init=', [startAt endAt], 50);
        for i = startAt:endAt
            manualSkip = false;
            progressBar(i);
            Did = Dids(i);

    %         Gid = dbLookup('Gid',  'Did', Did);
    %         plotSyncsAndSpikes(Gid(1));

            fprintf( [outOf(i, length(Dids)), '[Did = %d]'], Did );
            [nFramesPerPres_Runs, endStatus] = dbParseSyncs(Did);

            if (endStatus < 0)
                fprintf('Couldn''t parse syncs\n');
                continue;
            end
            if endStatus > 1
                3;
            end
            if strcmpi(stimType, 'Grating')
                [nDispDB, nSustDB, nPre, nPost, presNum, nFadeIn, nFadeOut] = getFieldsFromDatabaseTable(hnd, [frameFields, extraGratingFrameFields], stimTable, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
            else
                [nDispDB, nSustDB, nPre, nPost, presNum] = getFieldsFromDatabaseTable(hnd, frameFields, stimTable, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
                nFadeIn = 0;
                nFadeOut = 0;
            end
            nFadeIn = unique(nFadeIn);      assert(length(nFadeIn) == 1);
            nFadeOut = unique(nFadeOut);    assert(length(nFadeOut) == 1);
            nPre = unique(nPre);            assert(length(nPre) == 1);
            nPost = unique(nPost);          assert(length(nPost) == 1);
%             nSustDB = nSustDB);      assert(length(nSustDB) == 1);        
            nSustDB_av = median(nSustDB);
            nSustSyncs_av = median(nFramesPerPres_Runs)-nFadeIn-nFadeOut;
            if strcmpi(stimType, 'movie')
                nSustCorrect = dbGetNumberOfMoviePresFrames('Did', Did);
            elseif strcmpi(stimType, 'grating')
                nSustCorrect = ones(size(nSustDB))*roundToNearest( nSustSyncs_av, 10); % ***
                if nSustDB_av ~= nSustCorrect
                    keyboard;
                end            
            end
            assert( all( presNum == unique(presNum)) ); % make sure are unique, and in order.

%             nFramesFromSyncsCorrect = ones(size(nFramesPerPres_Runs)) * max(roundToNearest(nFramesPerPres_Runs, 1));
            nFramesFromSyncsCorrect = nFramesPerPres_Runs;
            nSustainedFrms_syncs = nFramesFromSyncsCorrect - nFadeIn - nFadeOut;
            nDisplayedFrms_syncs = nFramesPerPres_Runs + nPre + nPost;
            nDispCorrect = nDisplayedFrms_syncs; % ***

            indIncorrSust = nSustDB ~= nSustCorrect & ((~onlyCorrectZeros) | (nSustDB == 0));
            indIncorrDisp = nDispDB ~= nDispCorrect & ((~onlyCorrectZeros) | (nDispDB == 0));
            nIncorrSust = nnz( indIncorrSust );
            nIncorrDisp = nnz( indIncorrDisp );            
            
            if (nIncorrDisp == 0) && (nIncorrSust == 0)
                fprintf(' all correct \n');
            elseif (nIncorrDisp == 0) && (onlyCorrectnDisp && (nIncorrSust ~= 0));
                fprintf('[SKIPPING: %d/%d incorrect nsustained, but skipping.. ]] \n', ...
                    nIncorrSust, length(presNum) );                
            else
                if any(nDispDB > 0)
                    3; % [nDisplayedFrms_syncs  nDispDB]
                end
                indSustToFix = find(indIncorrSust);
                indDispToFix = find(indIncorrDisp);
                
                sustChanges = arrayfun(@(a,b) sprintf('x[%d->%d], ', a,b), nSustDB(indSustToFix), nSustCorrect(indSustToFix), 'un', false);
                [sustSum, nSust] = uniqueCount(sustChanges(:)');
                sustSummCell = [cellfun(@num2str, num2cell(nSust), 'un', 0); [sustSum(:)'] ];
                sustSummStr = [sustSummCell{:}];
                                
                dispChanges = arrayfun(@(a,b) sprintf('x[%d->%d], ', a,b), nDispDB(indDispToFix), nDispCorrect(indDispToFix), 'un', false);
                [dispSum, nDisp] = uniqueCount(dispChanges);
                dispSummCell = [cellfun(@num2str, num2cell(nDisp), 'un', 0); [dispSum(:)'] ];
                dispSummStr = [dispSummCell{:}];                

                dispStr = iff(~isempty(indDispToFix), ['nDisplayed: ' dispSummStr], '');
                sustStr = iff(~isempty(indSustToFix) && ~onlyCorrectnDisp, ['nSustained: ' sustSummStr], '');
                fprintf('Fixing %s %s\n', sustStr, dispStr)
                3;
                if manualSkip 
                    continue;
                end
                for pres_i = 1:length(presNum)
                    % make sure nSustained are correct
                    
                    if ~onlyCorrectnDisp 
                        if nSustDB(pres_i) ~= nSustCorrect(pres_i)
                            updateValueInDatabaseTable(hnd, nSustCorrect(pres_i), 'LNG_N_SUSTAINED_FRM', stimTable, {'DATAFILE_ID', Did; 'LNG_PRESENT_NO', pres_i});
                        end
                    end

                    if nDispDB(pres_i) ~= nDispCorrect(pres_i)
                        updateValueInDatabaseTable(hnd, nDispCorrect(pres_i), 'LNG_N_DISPLAYED_FRM', stimTable, {'DATAFILE_ID', Did; 'LNG_PRESENT_NO', pres_i});
                    end

                end
            end
            % make sure nDisplayed are correct

            3;
        end % of Dids for this stimulus
    
    end % of loop for all stimuli.
    
    
end
    
    

    
    
    

% fixSingleGratings = false;
% if fixSingleGratings
% 
%     load('cellsGroups_grating_dSng.mat');
%     singleGratDids = [gratingGroups_dSng(:).Did];
%     
%     [nDisp, nSust, nPre, nPost, presNum, Did, nFadeIn, nFadeOut] = getFieldsFromDatabaseTable(hnd, [frameFields, extraGratingFrameFields], stimTable, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');            
%     for i = 1:length(Dids)
%         Did = Dids(i);
%         Gid = dbLookup('Gid',  'Did', Did);
%         nSyncs = dbLookupNumSyncs(Did);
%         if ~mod(nSyncs, 10) 
%             keyboard;
%         end
% 
%         nSustI = nSust(i);
%         nDispI = nDisp(i);
%         nFadeInI = nFadeIn(i);
%         nFadeOutI = nFadeOut(i);
% 
%         if nSyncs > 0
%             nDisplayed_syncs = (nSyncs-1);
%             nSustained_syncs = nDisplayed_syncs - nFadeInI+nFadeOutI;
%     %         try
%             plotSyncsAndSpikes(Gid(1));
%     %         end
%             3;
% 
%             if ~(nSustained_syncs == nSustI) 
%                 updateValueInDatabaseTable(hnd, nSustained_syncs, 'LNG_N_SUSTAINED_FRM', 'TBL_GRATING_PRES', {'DATAFILE_ID', Did});
%             end
% 
%             if ~(nDisplayed_syncs == nDispI)                         
%                 updateValueInDatabaseTable(hnd, nDisplayed_syncs, 'LNG_N_DISPLAYED_FRM', 'TBL_GRATING_PRES', {'DATAFILE_ID', Did});                        
%             end
%         end
% 
%     end
% 

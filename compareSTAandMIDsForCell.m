function compareSTAandMIDsForCell(Gid, cellId)

    windows = {'best', [29, 62], [58, 91]};
    nWindows = length(windows);
    s = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
    
    showOddEven = 1;
%     M = iff(showOddEven, 4, 2);
    
    figure(1); clf;
    for i = 1:length(windows)
        if strcmp(windows{i}, 'best')
            STA{i} = s.STAs.STA;                 %#ok<*AGROW>
        else
            idx = find(floor(s.STAs.fixedWindowSTAs.windows_ms(:,1)) == windows{i}(1), 1);
            STA{i} = s.STAs.fixedWindowSTAs.STAs{idx};
        end        

        if ~showOddEven
            subplot(2,nWindows,i); 
            imagesc(mean(STA{i}, 3));
        else            
            STA_odd{i} = STA{i}(:,:,1);
            STA_even{i} = STA{i}(:,:,2);            
            subplot(4,nWindows,i);          imagesc(STA_odd{i});
            subplot(4,nWindows,i+nWindows); imagesc(STA_even{i});            
        end
                
    end
    
    3;
    %%
    for i = 1:length(windows);

        if showOddEven
            S_odd = load(mid_getPreferredMIDfile(Gid, cellId, windows{i}, 'odd'));   MID_odd{i} = S_odd.MID;
            S_even = load(mid_getPreferredMIDfile(Gid, cellId, windows{i}, 'even')); MID_even{i} = S_even.MID;
            subplot(4,nWindows,i+nWindows*2);  imagesc(MID_odd{i});
            subplot(4,nWindows,i+nWindows*3);  imagesc(MID_even{i});            
        else
            S = load(mid_getPreferredMIDfile(Gid, cellId, windows{i}, 'all')); 
            MID{i} = S.MID;
            subplot(2,nWindows,i+nWindows);  imagesc(MID{i});
        end        
        
    end
    3;
    allMID = [MID_odd{1}(:), MID_odd{2}(:), MID_odd{3}(:), MID_even{1}(:), MID_even{2}(:), MID_even{3}(:)];
    allSTA = [STA_odd{1}(:), STA_odd{2}(:), STA_odd{3}(:), STA_even{1}(:), STA_even{2}(:), STA_even{3}(:)];



end
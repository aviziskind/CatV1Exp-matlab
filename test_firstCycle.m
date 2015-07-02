global keepFirstCycle
[sGids, sCellIds] = getAllGids('s');
sNSpk = zeros(size(sGids));

uGids = unique(sGids);
for i = 1:length(uGids)
    sd = siteDataFor('Gid', uGids(i), 1);
    allCellIds = sd.cellIds(sd.cellIds > 0);
    allNSpk = sd.nSpikes(sd.cellIds > 0);
    for j = 1:length(allCellIds)
       idx = find(sGids == uGids(i) & sCellIds' == allCellIds(j));
       assert(length(idx) == 1)
       sNSpk(idx) = allNSpk(j);
    end
    
end

new_ord = ord(sNSpk, 'descend');
3;

doOSP = 1;
doPSTH = 1;

if doOSP

    for i = [7, 16, 26, 29, 38, 50, 64, 76]
    %%
        Gid = sGids(new_ord(i));
        cellId = sCellIds(new_ord(i));
    %     Gid = sGids(i);
    %     cellId = sCellIds(i);
        xlab = 'Temporal Phase';
        ylab = 'Spatial Freq.';


        keepFirstCycle = 1;
        s_withFirstCycle = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
        figure(1);  clf; imageOSP(s_withFirstCycle.OSP, 'pref:ori');
    %     str = sprintf('(%d) Gid = %d. cellId = %d', i, Gid, cellId);
        str = {};
        title({str{:}, 'All Trials'});
        xlabel(xlab); ylabel(ylab, 'fontsize', 10);
    %     size(decompress(s_withFirstCycle.OSP.R_full))

        keepFirstCycle = 0;
        s_withoutFirstCycle = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
        figure(2);  clf; imageOSP(s_withoutFirstCycle.OSP, 'pref:ori');
        title({str{:}, 'All Trials except 1st'});    
        xlabel(xlab); ylabel(ylab, 'fontsize', 10);
    %     size(decompress(s_withoutFirstCycle.OSP.R_full))


        figure(3); clf;
        [f1, R1_rscl] = estimateSpkRateFactor(s_withFirstCycle.OSP.R);
        [f2, R2_rscl] = estimateSpkRateFactor(s_withoutFirstCycle.OSP.R);
        s_diff = s_withoutFirstCycle;
        s_diff.OSP.R = (s_withFirstCycle.OSP.R/f1 - s_withoutFirstCycle.OSP.R /f2)*f1;    
        imageOSP(s_diff.OSP, 'pref:ori')
        title({str{:}, 'First trial'});
        xlabel(xlab); ylabel(ylab, 'fontsize', 10);


        3;
    end

elseif doPSTH
    
    for i =  [  7, 16, 26, 29, 38, 50, 64, 76]
    %%
        Gid = sGids(new_ord(i));
        cellId = sCellIds(new_ord(i));
    %     Gid = sGids(i);
    %     cellId = sCellIds(i);
        xlab = 'Temporal Phase';
        ylab = 'Spatial Freq.';

        
        [nOri, nSpf, nPh, nTrials, isCphFlashed] = getGratingStimType(Gid);
        g_type = getGratingStimType(Gid);
        nCycles = g_type.nUniqueSeq;
        nRep = g_type.nSeqRep;
                
        keepFirstCycle = 1;
        s_withFirstCycle = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
        
        
        R = s_withFirstCycle.OSP.R;
        R_full = s_withFirstCycle.OSP.R_full;
        R_os = mean(R, 3);
        [~, i_max] = maxElement(R_os);
        %%
        R_full_trials = reshape(decompress(R_full), [nOri, nSpf, nPh, nCycles, nRep]);        
        R_ph_trials = mean(R_full_trials(i_max(1), i_max(2), :, :, :), 5);
        
        R_ph_trials = reshape(R_ph_trials, [nPh, nCycles]);
%         R_ph_trials = gaussSmooth(R_ph_trials, 0, 1, 1);
        
        R_ph_trials_norm = R_ph_trials/max(R_ph_trials(:));
        R_ph_trials_norm_offset = bsxfun(@plus, -R_ph_trials_norm*.8, 1:nCycles);
        
        ph = [1:nPh]*360/nPh;
        figure(66);  clf; plot(ph, R_ph_trials_norm_offset);
        
        
        set(gca, 'ytick', 1:nCycles, 'xtick', [0:90:360]);
        xlim([0 360]);
        ylim([0 nCycles+.5])
        xlabel('Phase'); ylabel('Response per trial')
        axis ij
        %%
        
    %     str = sprintf('(%d) Gid = %d. cellId = %d', i, Gid, cellId);
%         str = {};
%         title({str{:}, 'All Trials'});
%         xlabel(xlab); ylabel(ylab, 'fontsize', 10);
    %     size(decompress(s_withFirstCycle.OSP.R_full))




        3;
    end
    
    
end

function plotSpikesWithStimulus(Gid, cellId)

    plotOri = 1;
    plotSpf = 1;
    showCorticalState = 0;
    
    spks = getSpikes(Gid, cellId, 'sec');
    nPtsPerStim = 2;
    
    [frameStimIds_ori, uOri, uSpf] = getStimulusFrameSequence(Gid, 'O');    
    [frameStimIds_spf]             = getStimulusFrameSequence(Gid, 'S');    
    
%     ori_color = hsv(length(uOri));
    
    sd = siteDataFor(Gid);
    nPres = length(sd.presIds);
    stimType = getGratingStimType(Gid);
    
    nStimTot = length(frameStimIds_ori);
    nStimsPerPres = nStimTot/nPres;
    
    
    syncTimes_sec0 = dbGetSyncs('Gid', Gid, 'sec');    
    dfInfo = sd.dataFileInfo;
    [nOri, nSpf, nPh] = deal(length(sd.ori_deg), length(sd.spPeriod_pix), length(sd.spPh_deg));
    nStim = nOri*nSpf*nPh;
            
    lastTick = dfInfo.duration_sec;

    syncTimes_sec = [0; syncTimes_sec0; lastTick]; % add my own ticks at start & end of experiment.
%     syncFrmIds = zeros(1,length(syncTimes_ticks));
    
    t_end = syncTimes_sec(end);

    if flashedOrDrifting(Gid) == 1
        frmLength = sd.frameLength_ms/1000;
    else
        frmLength = sd.tempPeriod_sec;
    end
    
    if plotOri || plotSpf
        cumulativeFrames = [0, [1:nPres]*nStimsPerPres];    
        frameIdForIntervalId  = zeros(1, length(syncTimes_sec)   );  

        for iPres = 1:nPres
            offset = iPres;
            frame_idxs = cumulativeFrames(iPres)+1:cumulativeFrames(iPres+1);
            frameIdForIntervalId(frame_idxs+offset) = frame_idxs;        
        end    
    
        pixEdges = [0:frmLength/nPtsPerStim:t_end];
        pixCenters = binEdge2cent(pixEdges);
        nPix = length(pixCenters);

        pixSyncIds = binarySearch(syncTimes_sec, pixCenters, [], -1);
        pixFrameIds = frameIdForIntervalId (pixSyncIds);
        idxn0 = pixFrameIds > 0;

        pixStimIds_ori = nan(1, nPix);    
        pixStimIds_ori(idxn0) = frameStimIds_ori(pixFrameIds(idxn0));    

        pixStimIds_spf = nan(1, nPix);    
        pixStimIds_spf(idxn0) = frameStimIds_spf(pixFrameIds(idxn0));    

        pixStimIds_ori_plot = pixStimIds_ori / max(pixStimIds_ori); pixStimIds_ori_plot(1) ;
        pixStimIds_spf_plot = pixStimIds_spf / max(pixStimIds_spf); pixStimIds_spf_plot(1) ;
        
    end

    nSpc = 5;    
    
    figure(1); clf;
    plot_i = 1;
    [hax_spk, hax_ori, hax_spf] = deal([]);
    
    if plotOri
        subplot(nSpc,1,plot_i);    
        plot_i = plot_i+1;
        oriLims = lims(pixStimIds_ori_plot);
        imagesc(pixCenters, 1,  pixStimIds_ori_plot);
        xlim([0, t_end]);
%         caxis([-1, length(uOri)]);
        caxis([0, 1]);
        hax_ori = gca;
        set(gca, 'ytick', [], 'xtick', []);
        ylabel('ori');
    end
    
    if plotSpf
        subplot(nSpc,1,plot_i);    
        plot_i = plot_i+1;
        imagesc(pixCenters, 1, pixStimIds_spf_plot);        
        xlim([0, t_end]);    
%         caxis([-1, length(uSpf)*2]);
        caxis([0, 1]);
        hax_spf = gca;
        set(gca, 'ytick', [], 'xtick', []);
        ylabel('spf');
    end

    if plotOri || plotSpf
        nCol = min(max(length(uOri), length(uSpf))*2, 255);
        colormap([0 0 0; hsv(nCol)]);
    end
    
    subplot(nSpc,1,plot_i:nSpc);
    binE = linspace(0, t_end, 300);    
    binC = binEdge2cent(binE);
    binVals = histcnt(spks, binE);
    bar(binC, binVals, 1);
    hax_spk = gca;
    xlim(binE([1, end]));
    
    h_ax = [hax_ori, hax_spf hax_spk];
    h_link = linkprop(h_ax, 'xlim');
    
    ylabel(sprintf('Group %d; cell %d', Gid, cellId))
    
    %%
    [recur_t, recur_nFrm] = getStimulusRecurrenceTime(Gid);
    [decor_t, decor_nFrm] = getStimulusDecorrelationTime(Gid);
    hLine_recur = drawVerticalLine(recur_t*[1:5], 'color', 'r', 'linewidth', 2, 'linestyle', '-');
    hLine_decor = drawVerticalLine(decor_t*[1:5], 'color', 'b', 'linewidth', 2, 'linestyle', ':');
    str_recur = sprintf(' Recurrence time = %.3f sec (%.1f stimuli)\n', recur_t, recur_nFrm);
    str_decor = sprintf(' Decorrelation time = %.3f sec (%.1f stimuli)\n', decor_t, decor_nFrm);
    fprintf(str_recur);
    fprintf(str_decor);
%     title(str)
%     text(mean(xlim), mean(ylim), str)
    if recur_t > binE(end) 
        xlim([binE(1), recur_t*1.01]);
    end
    legend([hLine_recur(1), hLine_decor(1)], {str_recur, str_decor}, 'location', 'best');
    
    if showCorticalState
        state_fun = getCorticalState(Gid, cellId);
        state_val = feval(state_fun, binC);
        hold on;
        plot(binC, state_val, 'r');
        
    end
    
    
    3;


end
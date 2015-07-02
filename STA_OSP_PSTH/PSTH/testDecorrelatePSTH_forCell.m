function testDecorrelatePSTH_forCell(Gid, cellId)


%     redo = true;
%     if ~exist('decorrelationTestStart.mat', 'file') || redo
%     [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms);
%     nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;

    frameLength_ms = getFrameLength('Gid', Gid);
    OSP_30_60 = getOriSpfPhaseProfile_simple(Gid, cellId, [30 60], []);
    [PSTHdata, spkTsRelToFrame_ms] = getPSTHforCell(Gid, cellId, true); warning('dc:getPSTHformat', 'no longer correct format');
    timeWindow_ms = PSTHdata.timeWindow;
    OSP_best1PSTH = getOriSpfPhaseProfile_simple(Gid, cellId, timeWindow_ms, PSTHdata.windowProfile );        
    [nOri, nSp, nPh] = size(OSP_30_60); %#ok<NASGU>
                        
%         save('decorrelationTestStart.mat', '*');
%     else
%         load('decorrelationTestStart.mat');
%     end

  
    [dcOSP, corOSP] = getOriSpfPhaseProfile_decor(Gid, spkTsRelToFrame_ms, timeWindow_ms, frameLength_ms);

    if isJavaRunning
        fig_id = 167;
        figlabel = ['Group ' num2str(Gid) '.  Cell ' num2str(cellId)];
        figname =  ['Group' num2str(Gid) '_Cell' num2str(cellId)];
        
        figure(fig_id);
        subplot(1,4,1); imageOSP(OSP_30_60);        title('30 - 60');
        ylabel(figlabel);
        subplot(1,4,2); imageOSP(OSP_best1PSTH);    title( sprintf('%d - %d, ', round(timeWindow_ms(1)), round(timeWindow_ms(2)) ) );
        subplot(1,4,3); imageOSP(corOSP);           title('full PSTH (correlated)');
        subplot(1,4,4); imageOSP(dcOSP);            title('decorrelated');
        saveas(fig_id, [figname '.fig'])
        
    end
    
    filename = ['Group' num2str(Gid) '_Cell' num2str(cellId) '.mat'];
    save(filename);
    
end
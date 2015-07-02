function plotPhaseTCcc_vs_F1oDC

    [fd, gratingType_s] = curGratingType;
    S = load([gratingType_s 'GratingCells_all.mat']);
    allCells = S.allCells;
    nCells = length(allCells);
    allStats = [allCells.stats];
    isRep = [allStats.isRep];
    
    allStimPhaseTcCCs = zeros(1, nCells);
    allStimPhaseTcDots = zeros(1, nCells);
    F1oDCs = zeros(1,nCells);
    
    progressBar('init-', nCells, 40);
    for i = 1:nCells
%                 fprintf('*');
        progressBar;                
        Gid = allCells(i).Gid;
        cellId = allCells(i).cellId;
        th = 0.80;
                
%                 [tgtWind_osp, phaseTC_CCs, stimF1oDCs] = getOspDataForPsthWindow(Gid, cellId, [], [], L, R, {'osp', 'phaseTC_CCs', 'stimF1oDCs'});
        if fd == 1
            [L_bin, R_bin] = dealV(allCells(i).PSTH.timeWindow_bins);
            windowProfile = allCells(i).PSTH.windowProfile;
        else
            [L_bin, R_bin] = deal(1);
            windowProfile = [];
        end
        [tgtWind_osp, phaseTC_CCs, phaseTC_CC_dots, stimF1oDCs] = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, windowProfile, {'osp', 'phaseTC_CCs', 'phaseTC_Dots', 'stimF1oDCs'});

        % x_axis: mean f1odc of stim > 80%                            
        F1oDCs(i) = getAllAboveTh(stimF1oDCs, tgtWind_osp, th);                
                
        % y axis: cc of phase tuning cc's & windows.                
        allStimPhaseTcCCs(i) = getAllAboveTh(phaseTC_CCs, tgtWind_osp, th);                
        allStimPhaseTcDots(i) = getAllAboveTh(phaseTC_CC_dots, tgtWind_osp, th);                
                
%                 tgtWind_osp_norm = tgtWind_osp / max(tgtWind_osp(:));
%                 ph_ccs_wgt = phaseTC_CCs .* tgtWind_osp_norm;                
%                 allStimPhaseTcCCs(i) = doPearsonCorr(tgtWind_osp_norm(:), ph_ccs_wgt(:));
    end
%     return;
    xLab = sprintf('F1/DC (averaged over stimuli %d%%+ of max)', th*100);
    
    figure(243); clf; hold on; box on;
    plot(F1oDCs(isRep), allStimPhaseTcCCs(isRep), 'bo');
    plot(F1oDCs(~isRep), allStimPhaseTcCCs(~isRep), 'ro');
    title('CC');
    xlabel(xLab);    
    ylabel( sprintf('phase-TC CC (av. over stimuli %d%%+ of max)', th*100) );
    xlim([0 2]); ylim([-1 1]); 

    figure(244); clf; hold on; box on;
    plot(F1oDCs(isRep), allStimPhaseTcDots(isRep), 'bo');
    plot(F1oDCs(~isRep), allStimPhaseTcDots(~isRep), 'ro');
    title('Dot');
    xlabel(xLab);    
    ylabel( sprintf('phase-TC Dots (av. over stimuli %d%%+ of max)', th*100) );
    xlim([0 2]); ylim([0 1]); 
end

function y = getAllAboveTh(X, osp, th)
    idx = (osp(:) > max(osp(:))*th);
    y = mean(X(idx));
end

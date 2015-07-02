function tmp_compareNStdFlashedDrifting

    doHistNStd = 0;
    doHistBckg = 0;
    doHistR_max = 1;
    
    if doHistNStd
        figure(11);
        S = load('driftingGratingCells_GLFcuw8_degree_all.mat');
        allStats = [S.allCells.stats];
        isOri = arrayfun(@(s) isfield(s.tuningStats, 'oriStats_si'), allStats);
        allOriTunStats = [allStats(isOri).tuningStats];
        allOriStats_si = [allOriTunStats.oriStats_si];
        nStds = [allOriStats_si.nStdAboveBckg];
        nTot = length(nStds);
        nStds(nStds > 175) = 200;
        binE = [-1.5:1.5:201];
        binC = binEdge2cent(binE);
        binV = histcnt(nStds, binE);
        bar(binC, binV, 1)
        drawVerticalLine(3, 'color', 'r')
        xlim([-5, 202])
        nAbove3 = nnz(nStds > 3);
        xlabel('N Std deviations above bckg')
        length(nStds)
        title(sprintf('Drifting Gratings (N = %d). %.1f%% above 3', nTot, nAbove3/nTot*100));


        figure(12);
        S = load('flashedGratingCells_GLFcuw8_degree_all');
        allStats = [S.allCells.stats];
        allTunStats = [allStats.tuningStats];
        allOriStats_si = [allTunStats.oriStats_si];
        nStds = [allOriStats_si.nStdAboveBckg];

        nStds(nStds > 15) = 15;
        nTot = length(nStds);
        binE = [-1:.25:15];
        binC = binEdge2cent(binE);
        binV = histcnt(nStds, binE);
        bar(binC, binV, 1)
        xlim([-1, 15])
        drawVerticalLine(3, 'color', 'r')
        nAbove3 = nnz(nStds>3);
        xlabel('N Std deviations above bckg')
        title(sprintf('Flashed Gratings (N = %d). %.1f%% above 3', nTot, nAbove3/nTot*100));
    end
    
    
    if doHistBckg        
        
        [dgGids, dgCellIds] = getAllGids('o'); nD = length(dgGids);
        [bckgMeans_d, bckgStds_d, meanRates_d] = deal(zeros(1, nD));
        
        
        progressBar('init-', nD)
        for i = 1:nD
            [~,~, meanRates_d(i), bckgSamples] = dbGetCellSpkStimHists(dgGids(i), dgCellIds(i));
            bckgMeans_d(i) = mean(bckgSamples);
            bckgStds_d(i) = std(bckgSamples);
%             bckgSamples_d{i}
            progressBar(i);
        end            
        bckgMeans_d_orig = bckgMeans_d;
        
        bckgMax = 15;
        bckgMeans_d(bckgMeans_d > bckgMax) = bckgMax;
        figure(13); hist(bckgMeans_d, 41);
        title(sprintf('Drifting Gratings (N = %d). Mean = %.2f', nD, mean(bckgMeans_d))); 
        xlabel('Background firing rate (Hz)')

        figure(21); hist(meanRates_d, 41);
        title(sprintf('Drifting Gratings (N = %d). Mean = %.2f', nD, mean(meanRates_d))); 
        xlabel('Mean stimulus-induced firing rate (Hz)')
                
        
        [fgGids, fgCellIds] = getAllGids('f'); nF = length(fgGids);     
        [bckgMeans_f, bckgStds_f, meanRates_f] = deal(zeros(1, nF));
        progressBar('init-', nF)
        for i = 1:nF
            [~,~, meanRates_f(i), bckgSamples] = dbGetCellSpkStimHists(fgGids(i), fgCellIds(i));
            bckgMeans_f(i) = mean(bckgSamples);
            bckgStds_f(i) = std(bckgSamples);
            progressBar(i);
        end            
        bckgMeans_f_orig = bckgMeans_f;
        
        figure(14);
        bckgMeans_f(bckgMeans_f > bckgMax) = bckgMax;        
        hist(bckgMeans_f, 41);                
        title(sprintf('Flashed Gratings (N = %d). Mean = %.2f', nF, mean(bckgMeans_f) )); xlabel('Background firing rate (Hz)')
        
        figure(22); hist(meanRates_f, 41);
        title(sprintf('Flashed Gratings (N = %d). Mean = %.2f', nF, mean(meanRates_f))); 
        xlabel('Mean stimulus-induced firing rate (Hz)')
    3;
        
    end
    
    if doHistR_max
        S_d = load('driftingGratingCells_GLFcuw8_phase.mat');
        nD = length(S_d.allCells);
        R_maxs_d = zeros(1, nD);
        for i = 1:nD
            R = mean(S_d.allCells(i).R, 3);
            R_maxs_d(i) = max(R(:));
        end
        maxR = 100;
        R_maxs_d_orig = R_maxs_d;
        
        figure(15);
        R_maxs_d(R_maxs_d>maxR) = maxR;
        hist(R_maxs_d, 40);
        title(sprintf('Drifting Gratings (N = %d). Mean = %.2f', nD, mean(R_maxs_d_orig) )); xlabel('Peak firing rate (Hz)');
        
        S_f = load('flashedGratingCells_GLFcuw8_phase.mat');
        nF = length(S_f.allCells);
        R_maxs_f = zeros(1, nF);
        for i = 1:nF
            R = mean(S_f.allCells(i).R, 3);
            R_maxs_f(i) = max(R(:));
        end
        R_maxs_f_orig = R_maxs_f;
        figure(16);
        R_maxs_f(R_maxs_f>maxR) = maxR;
        hist(R_maxs_f, 40);
        title(sprintf('Flashed Gratings (N = %d). Mean = %.2f', nF, mean(R_maxs_f_orig) )); xlabel('Peak firing rate (Hz)');
        
        3;
    end
               
        
        
    
    
    
    
    
end
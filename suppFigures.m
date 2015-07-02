function suppFigures
doFigs = [4];

% Supplementary Figure 1c
% Gid = 652;
% figure(3);




% Supplementary Figure 4 : F1/DC ratios

    if any(doFigs == 4)
        for i = 1:2
            if i == 1
                title_str = 'Drifting Gratings';
                S = load('driftingGratingCells_GLFcuw8_degree.mat');
        
                stimTypes = {S.allCells.stimType};
                idx_spfBatch = strncmp(stimTypes, 'Grating:Spatial', 12);
            elseif i == 2
                title_str = 'Flashed Gratings';
                S = load('flashedGratingCells_GLFcuw8_degree.mat');
                        
                idx_spfBatch = 1:length(S.allCells);                
            end
                
            allF1oDCs = [S.allCells(idx_spfBatch).F1oDC_maxR_avP];
            figure(400+i);            
            idx_complex = allF1oDCs < 1;
            idx_simple = allF1oDCs >= 1;

            binE = linspace(0,2,19);
            binC = binEdge2cent(binE);
            binV_simp = histcnt(allF1oDCs(idx_simple), binE);
            binV_comp = histcnt(allF1oDCs(idx_complex), binE);        

            hb = bar(binC, [binV_comp(:), binV_simp(:)], 1, 'stacked');
            set(hb(2), 'facecolor', 'r');
            legend({'Complex', 'Simple'}, 'location', 'NE');
            xlabel('F1/DC'); ylabel('Number of Cells');
            title(title_str);

        end
    end


% Supplementary Figure 5 : w_ori vs DSI


end
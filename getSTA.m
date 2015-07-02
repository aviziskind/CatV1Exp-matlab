function getSTA(Gid, cellId, opt)

    persistent allSTAs
    timeWindow, trialMode, jackknifeMethod


    addFittedGaborStats = 1;
    addSTAs = 0;
        flipMIDforPosCorrWithSTA = 0;
        onlyIncludeMIDsWithFits = 0;
        
    addMID_fit = 0;
        addSigPixels_gabor = 1;
        addSigPixels_2stdMID = 0;
    
    addMID_jacks = 1;
    
    allTimeWindows = {'best', [29, 62], [58, 91]};
    allTrialModes = {'all', 'odd', 'even'};
    
    nWindows = length(allTimeWindows);
    nTrialModes = length(allTrialModes);
    
    
    S_cells = load('usedCells.mat');    
    getAllCells = 1;    
    if getAllCells         
        S.allGids = S_cells.allGids;
        S.allCellIds = S_cells.allCellIds;
    else
        S.allGids = S_cells.usedGids;
        S.allCellIds = S_cells.usedCellIds;    
    end
    
    nCells= length(S.allCellIds);
%     nCells = 4;
%     S = struct;
%     S.allMIDs = cell(nWindows,nCells);
%     S.allMIDs_odd = cell(nWindows,nCells);
%     S.allMIDs_even = cell(nWindows,nCells);
%     S.rsqr = zeros(1,nCells);
%     S.jackMeanCC = zeros(1,nCells);
%     S.gparams = cell(1,nCells);
%     S.MID_fit = cell(1,nCells);
%     S.above2std_mid = cell(1,nCells);
%     S.above2std_gabor = cell(1,nCells);
    
%     S.frac_above2std_mid = zeros(1,nCells);
%     S.frac_above2std_gabor = zeros(1,nCells);
    
    if addSTAs        
        S_indivCells = load( getFileName('indiv', 'movie_fg') );         
        S.allSTAs = cell(1,nCells);        
    end
    
    sig_pixels_gabor_flds = {}; %#ok<NASGU>
    if addSigPixels_gabor
        sig_pixels_gabor_flds = { 'above2std_gabor', [], 'frac_above2std_gabor', [] };
    end
    
    sig_pixels_mid_flds = { };
    if addSigPixels_2stdMID 
        sig_pixels_mid_flds = { 'above2std_mid', [], 'frac_above2std_mid', [] };        
    end
    
    mid_jacks_flds = {}; %#ok<NASGU>
    if addMID_jacks
        mid_jacks_flds = {'MID_jacks', []};
    end
            
    
    emptyStruct = struct('MID', [], 'jackMeanCC', [], 'rsqr', [], 'gparams', [], 'MID_fit', [], sig_pixels_gabor_flds{:}, sig_pixels_mid_flds{:}, mid_jacks_flds{:}, 't_calc', []);
    







end




function calculateSTA(Gid, cellId, timeWindow, trialMode, jackknifeMethod) 

        if any(strcmp('STA', dataToCalc))                                    
            r_ph = decompress (s.osp_ph{w_idx});                         
            STA = getSTAfromOSP(Gid, r_ph);            
            s.STA{w_idx} = single(STA);                        
        end
        
        if any(strcmp('STA_oe', dataToCalc))
            r_ph_oe = decompress (s.osp_ph_oe{w_idx});
            STA_odd = getSTAfromOSP(Gid, r_ph_oe(:,:,:,1));
            STA_even = getSTAfromOSP(Gid, r_ph_oe(:,:,:,2));
            s.STA_oe{w_idx} = cat(3, single(STA_odd), single(STA_even));                        
        end
        
        stimInfo = getGratingStimType(Gid);
        nJackSegments_STA = switchh(stimInfo.nTrials, [4, 16,   10], [4, 4,  5]);
        jackknifeMethod_STA = 'trials';
        
        if any(strcmp('STA_jack', dataToCalc))
            r_full = decompress (s.osp_full{w_idx});
            s.STA_jack{w_idx} = getJackknifedSTA(r_full, Gid, 'all', nJackSegments_STA, jackknifeMethod_STA);
        end

        if any(strcmp('STA_odd_jack', dataToCalc))
            r_full = decompress (s.osp_full{w_idx});
            s.STA_odd_jack{w_idx} = getJackknifedSTA(r_full, Gid, 'odd', nJackSegments_STA, jackknifeMethod_STA);
        end

        if any(strcmp('STA_even_jack', dataToCalc))
            r_full = decompress (s.osp_full{w_idx});
            s.STA_even_jack{w_idx} = getJackknifedSTA(r_full, Gid, 'odd', nJackSegments_STA, jackknifeMethod_STA);
        end
        
        
        
        
        
        
end
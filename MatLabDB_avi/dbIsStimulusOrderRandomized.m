function isRandOrder = dbIsStimulusOrderRandomized(Gid)
%     S_amend = dbGetSyncAmendments(idType, idVal)        
%     S_amend = dbGetSyncAmendments(Gid)

    persistent allExpStimulusOrder saveCount 
    
    stimOrder_file = [CatV1Path 'MatLabDB_avi' filesep 'allExpStimulusOrder.mat'];
        
    redo_all = 0;
    redo_current = 1;
    saveCountSpacing = 50;

    if strcmp(Gid, 'save')
        save(stimOrder_file, 'allExpStimulusOrder', '-v6');                
        saveCount = 0;
        return;
    end        

    
    if isempty(allExpStimulusOrder)        
        if exist(stimOrder_file, 'file') && ~redo_all
            S_file = load(stimOrder_file);
            allExpStimulusOrder = S_file.allExpStimulusOrder;
        else
            allExpStimulusOrder = struct;
        end        
        saveCount = 0;        
    end
        
    grp_fld_name = sprintf('stimOrder_Gid_%d', Gid);        
        
    if (~isfield(allExpStimulusOrder, grp_fld_name) || redo_current)
        [isRandOrder, L] = calcIfStimulusOrderRandomized(Gid);
        
        allExpStimulusOrder.(grp_fld_name) = L;
        saveCount = saveCount + 1;

        if saveCount > saveCountSpacing
            save(stimOrder_file, 'allExpStimulusOrder', '-v6');        
            saveCount = 0;
        end                
    end
    
%     if nargout > 0
        isRandOrder = allExpStimulusOrder.(grp_fld_name);
%     end


end

function [isRandOrder, L] = calcIfStimulusOrderRandomized(Gid)
    
    sd = siteDataFor(Gid);
    if flashedOrDrifting(Gid) == 1
        isRandOrder = 1;
        return;
    end
    
    Did = sd.Did;
    presIds = sd.presIds;
    nPres = length(presIds);
    allOris = sd.ori_deg;
    allSpfs = sd.spPeriod_pix;
    
%     getStimulusFrameSequence
    
    presOri = zeros(1,nPres);
    presSpf = zeros(1,nPres);
    hnd = dbOpenExpDb;
    for pres_i = 1:nPres
        [phase_deg, sp_pix, ori_deg] = getGratingOSPs_grating(hnd, presIds(pres_i));
        presOri(pres_i) = ori_deg(1);
        presSpf(pres_i) = sp_pix(1);
    end
    presOri_idx = binarySearch(allOris, presOri);
    presSpf_idx = binarySearch(allSpfs, presSpf);
    
    assert(all(presOri_idx ~= 0))
    assert(all(presSpf_idx ~= 0))
    
    smoothFactors = [0, 2, 3, 5, 10]; nS = length(smoothFactors);
    L = cell(1,nS+1);
    for sm_i = 1:nS
        diffPresOri_sm = gaussSmooth( presOri_idx-mean(presOri_idx), smoothFactors(sm_i), [], 1);
        diffPresSpf_sm = gaussSmooth( presSpf_idx-mean(presSpf_idx), smoothFactors(sm_i), [], 1);
    
        L{sm_i}(1) = mean( abs(diffPresOri_sm) );
        L{sm_i}(2) = mean( abs(diffPresSpf_sm) );
    end
%      diff(presOri_idx)
    
    
    diff(diffPresOri_sm)
    idx1 = 1:nPres-1;
    idx2 = 2:nPres;
    
    ori_avD_all = circDist(presOri_idx(idx1), presOri_idx(idx2), length(allOris));
    spf_avD_all = circDist(presSpf_idx(idx1), presSpf_idx(idx2), length(allSpfs));

    L{nS+1}(1) = mean(ori_avD_all);
    L{nS+1}(2) = mean(spf_avD_all);
    
      
%     ori_th = 1.3; % unordered is up to ~1.12.  have ~3 ordered at 1.5, and most ordered are 
%     ori_smoothFactor = 10; % unordered is up to ~1.12.  have some 
    
    
    
    
    
end


%{
allGids = [getAllGids('d')];
L = cell(1,length(allGids));
progressBar('init-', length(allGids), 60);
for i = 1:length(allGids)
    L{i} = dbIsStimulusOrderRandomized(allGids(i));
    progressBar(i);
end

%}

%{
L3 = cellfun(@(C) C{3}, L, 'un', 0);
L_cat3 = cat(1, L3{:});

L5 = cellfun(@(C) C{4}, L, 'un', 0);
L_cat5 = cat(1, L5{:});

L10 = cellfun(@(C) C{5}, L, 'un', 0);
L_cat10 = cat(1, L10{:});

%}
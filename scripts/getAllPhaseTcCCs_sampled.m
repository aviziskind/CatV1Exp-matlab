function [ph_ccs, ph_cc_ps] = getAllPhaseTcCCs_sampled(osp_trials, corrType, trialDivisionType)
    nMax = 1;

    if nargin < 2
        corrType = 'cc';
    end
    CC = 1; DOT = 2; 
    corrType = switchh(corrType, {'cc', 'dot'}, {CC, DOT});

    
    if nargin < 3
        trialDivisionType = 'odd_vs_even';
    end
    
%     justDo_odd_vs_even = 1;
    
    doTests = false;
    scrambleTrialsForMorePerms = false;
    doAllPvalsAtEnd = true;
    
    [nOri, nSp, nPh, nTrials] = size(osp_trials);    
    
    if strcmp(trialDivisionType, 'odd_vs_even')
        halfTrials1 = odd([1:nTrials]);
        halfTrials2 = ~halfTrials1;
        
        perms = [find(halfTrials1), find(halfTrials2)];
%         perms = [1:2:nTrials, 2:2:nTrials];

        nPerms = 1;
        
    elseif strcmp(trialDivisionType, 'half_odd_vs_half_even')
    
%         perms = [1:2:nTrials, 2:2:nTrials];
        halfTrials1 = odd(floor([1:nTrials] / 2));
        halfTrials2 = ~halfTrials1;
               
        % note: would have been "neater" to use 'ceil' instead of 'floor' for
        % half-odd/even definition: 
        %  ceil:  half-odd = [1,2,  5,6]. half-even = [3,4,  7,8]
        %  floor: half-odd = [1, 4,5, 8,9]. half-even = [2,3,  6,7]
        % but i'm using floor here for the case where there are 10 trials.
        % so that both 'half-odd' and 'half-even' have 5 trials 
        % (with ceil, 'half odd' has 6 trials, 'half-odd' has 4
        % trials).
        if ~odd(nTrials)
            assert( nnz(halfTrials1) == nnz(halfTrials2));
        else
            assert( abs(nnz(halfTrials1) - nnz(halfTrials2)) == 1);
        end
        
        perms = [find(halfTrials1), find(halfTrials2)];
%         perms = [1:2:nTrials, 2:2:nTrials];

        nPerms = 1;
    elseif strcmp(trialDivisionType, 'first_vs_second_half')
    
%         perms = [1:2:nTrials, 2:2:nTrials];
        halfTrials1 = [1:nTrials] <= nTrials/2;
        halfTrials2 = ~halfTrials1;
        
        perms = [find(halfTrials1), find(halfTrials2)];
%         perms = [1:2:nTrials, 2:2:nTrials];

        nPerms = 1;        
    else
        if nTrials <= 16  % for flashed grating experiments
            perms = nchoosek(1:nTrials, nTrials/2);    
            nPerms = size(perms,1)/2;

            perms = [perms(1:nPerms,:), perms(nPerms*2:-1:nPerms+1,:)];
                p = randperm(nPerms);        
                perms = perms(p, :);        

        else              % drifting grating experiments        
            nPerms = nMax;
            [tmp, perms] = sort(rand(nPerms, nTrials), 2);        
        end
    end        

    if doTests
        assert( all (sum(perms,2)  == sum(1:nTrials) ) );
    end
    
    scrambleTrials = false;    
    nPermsToDo = nMax;
    
    if nPerms <= nMax
        if scrambleTrialsForMorePerms % can shuffle around for each individual phase        
            nPermsWithPh = nPerms^nPh;
            nPermsToDo = min(nPermsToDo, nPermsWithPh);            

            ph_perm_idx = randi(nPerms, nPermsToDo, nPh);
            scrambleTrials = true;
        else
            nPermsToDo = nPerms;    
        end        
    end
        

%     n_total = nchoosek(nTrials, nTrials/2);
    
%     stimPSTH_vals_allTrials = reshape(stimPSTH_vals_allTrials, [nBins, nOri, nSp, nPh, nTrials])
%     osp_trials = sum(osp_trials([L_bin_wind_oe : R_bin_wind_oe], :, :, :, :), 1);
    [ccs_smp, ph_cc_ps_smp] = deal( zeros(nOri,nSp, nPermsToDo) ); 
    calcPvals = nargout > 1;
    idx1 = 1:nTrials/2; idx2 = setdiff(1:nTrials, idx1); %nTrials/2+1:nTrials;
    
    
    if doTests
        nTot_ph = squeeze(sum(osp_trials,4));
        nTot    = sum( nTot_ph, 3);
    end    
    
    
    ph_tc1 = zeros(nPh, 1);
    ph_tc2 = zeros(nPh, 1);
%%
    for prm_i = 1:nPermsToDo
%         trialOrder = randperm(nTrials);
%         trialOrder1 = trialOrder(1:nTrials/2);
%         trialOrder2 = trialOrder(nTrials/2+1:nTrials);
        if ~scrambleTrials
%             trialOrder1 = idx1;
%             trialOrder2 = idx2;
            trialOrder1 = perms(prm_i, idx1);
            trialOrder2 = perms(prm_i, idx2);
        else 
%             trialOrder1 = cell(1,nPh);
%             trialOrder2 = cell(1,nPh);
%             for ph_i = 1:nPh
%                 trialOrder1{ph_i} = perms(ph_perm_idx(prm_i,ph_i), idx1);
%                 trialOrder2{ph_i} = perms(ph_perm_idx(prm_i,ph_i), idx2);                
%             end  
3;
            trialOrder1 = arrayfun(@(ph_i) perms(ph_perm_idx(prm_i,ph_i), idx1), 1:nPh, 'un', 0);
            trialOrder2 = arrayfun(@(ph_i) perms(ph_perm_idx(prm_i,ph_i), idx2), 1:nPh, 'un', 0);
        end
        
        for oi = 1:nOri
            for si = 1:nSp
                if ~scrambleTrials
                    ph_tc1 = sum(osp_trials(oi, si, :, trialOrder1), 4);
                    ph_tc2 = sum(osp_trials(oi, si, :, trialOrder2), 4);
                else
                    for ph_i = 1:nPh
                        ph_tc1(ph_i) = sum(osp_trials(oi, si, ph_i, trialOrder1{ph_i}), 4);
                        ph_tc2(ph_i) = sum(osp_trials(oi, si, ph_i, trialOrder2{ph_i}), 4);
                    end
                end
                
                if doTests
                    tot_ph = round(ph_tc1(:) + ph_tc2(:)) ;
                    tot = round(sum(tot_ph)) ;
                    
                    Tot_ph = round( squeeze(nTot_ph(oi, si,:)) );
                    Tot    = round(nTot(oi, si));
                    assert( all( tot_ph == Tot_ph) );                    
                    assert( tot == Tot);                    
                end
                
                if corrType == CC
                    if (sum(ph_tc1) == 0) || (sum(ph_tc2) == 0)
                        ccs_smp(oi,si, prm_i) = nan;
                    else                    
                        if calcPvals && ~doAllPvalsAtEnd
                            [ccs_smp(oi,si, prm_i), ph_cc_ps_smp(oi,si,prm_i)] = ...
                                doPearsonCorr(ph_tc1, ph_tc2 );
                        else
                            [ccs_smp(oi,si, prm_i)] = ...
                                doPearsonCorr(ph_tc1, ph_tc2);
                        end
                    end
                elseif corrType == DOT
                    [ccs_smp(oi,si, prm_i)] = ...
                        normDotProd(ph_tc1, ph_tc2);                    
                end
            end
        end
    end
    ph_ccs = mean(ccs_smp, 3);
    if calcPvals && doAllPvalsAtEnd
        ph_cc_ps = pearsonPval('r', ph_ccs, nPh);
    end
            
end





%{


%%
N = 40;
trial_ids = 1:N;
tf_odd = odd(trial_ids);
tf_even = ~odd(trial_ids);
tf_new1 = odd(floor(trial_ids / 2));
tf_new2 = ~tf_new1;

A = [nnz(tf_odd & tf_new1), nnz(tf_even & tf_new1);
     nnz(tf_odd & tf_new2), nnz(tf_even & tf_new2)]

%}
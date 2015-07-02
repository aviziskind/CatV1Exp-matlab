function [ph_ccs, ph_cc_ps] = getAllPhaseTcCCs(R1, R2)
            
    [nOri, nSp, nPh] = size(R1);    
    assert(isequal(size(R1), size(R2)));
        
%     osp_trials = sum(osp_trials([L_bin_wind_oe : R_bin_wind_oe], :, :, :, :), 1);
    [ccs_smp, ph_cc_ps_smp] = deal( zeros(nOri,nSp) ); 
    calcPvals = nargout > 1;

    doAllPvalsAtEnd = true;

%%

    for oi = 1:nOri
        for si = 1:nSp
            ph_tc1 = squeeze(R1(oi, si, :));
            ph_tc2 = squeeze(R2(oi, si, :));

            if all(ph_tc1==0) || all(ph_tc2==0)
                ccs_smp(oi,si) = nan;
            else
                if calcPvals && ~doAllPvalsAtEnd
                    [ccs_smp(oi,si), ph_cc_ps_smp(oi,si)] = ...
                        doPearsonCorr(ph_tc1, ph_tc2 );
                else
                    [ccs_smp(oi,si)] = ...
                        doPearsonCorr(ph_tc1, ph_tc2);
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
function [stimF1oDCs, stimF1oDC_jackStds] = getStimF1oDCs(R_full)
    [nOri, nSpf, nPh, nTrials] = size(R_full);
    nStim = nOri * nSpf;
    calcJackStd = nargout > 1;
    if calcJackStd && nTrials <= 2
        error('Need all trials to calculate jackStdErrs')
    end
            
    r_stim = mean(reshape(R_full, nStim, nPh, nTrials), 3);
    
    stimF1oDCs = getF1oDC(r_stim);
    stimF1oDCs = reshape(stimF1oDCs, [nOri, nSpf]);
    
    if calcJackStd
        R_full_jackknifeTrials = getPhaseTuningJackknifedTrials(R_full, 'aa');
        stimF1oDC_jackknives = cellfun(@getF1oDC, R_full_jackknifeTrials);

        stimF1oDC_jackStds = zeros(nStim, 1);
        for si = 1:nStim        
            stimF1oDC_jackStds(si) = jackknifeStdErr(stimF1oDC_jackknives(si,:), stimF1oDCs(si));
        end        
        stimF1oDC_jackStds = reshape(stimF1oDC_jackStds, [nOri, nSpf]);
    end
            
end



%{
    nPh = size(R_full, 3); phs = linspace(0, 360, nPh+1); phs = phs(1:nPh);
    nTrials = size(R_full, 4);
    phaseTC_atPref_trials = squeeze( R_full(ori_pref_idx, spf_pref_idx, :, :) );
    assert(isequal(size(phaseTC_atPref_trials), [nPh, nTrials]));
    phaseTC_atPref = mean( phaseTC_atPref_trials, 2);
        
    F1oDC = getF1oDC(phs, phaseTC_atPref(:), 360);
    opt.F1oDC_calculated = F1oDC;
    
%     phaseTC_atPref_jacks = arrayfun(@(i) mean(  phaseTC_atPref_trials(:, setdiff(1:nTrials, i) )  ,2), 1:nTrials, 'un', 0);
    phaseTC_atPref_jacks = jackknifeAverageTrials( phaseTC_atPref_trials, 2 );
    
    
    %%
    if opt.doJackKnifeStdErr
        F1oDC_jacks = deal(  zeros(1,nTrials) );
        for jack_i = 1:nTrials
            F1oDC_jacks(jack_i) = getF1oDC(phs, phaseTC_atPref_jacks{jack_i}(:), 360);
        end
        F1oDC_stderr_jack = jackknifeStdErr(F1oDC_jacks, F1oDC);
        opt.F1oDC_stderr_calculated = F1oDC_stderr_jack;
    end
    
%}
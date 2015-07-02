function R_full_jackknifeTrials = getPhaseTuningJackknifedTrials(R_full, phaseOEmode)

            
    [n_ori, n_spf, n_ph, n_trials] = size(R_full);
    n_stim = n_ori*n_spf;
    R_full_stim_ph_trial = reshape(R_full, [n_stim, n_ph, n_trials]);
            
%             R_full_jackknifeTrials = cell(n_stim, n_trials);
%             for jj = 1:n_stim
%                 R_full_jackknifeTrials(jj,:) = jackknifeAverageTrials(R_full_stim_ph_trial(jj,:,:), 3);
%             end
%         
    if strcmp(phaseOEmode, 'aa')
        R_full_jackknifeTrials_stim = jackknifeAverageTrials(R_full_stim_ph_trial, 3);

        R_full_jackknifeTrials = cell(n_stim, n_trials);
        for jj = 1:n_stim
            R_full_jackknifeTrials(jj,:) = cellfun(@(x) x(jj,:), R_full_jackknifeTrials_stim, 'un', 0);
        end

        R_av_jack = mean(cat(3, R_full_jackknifeTrials_stim{:,:}),3);
        diffs = abs(R_av_jack - mean(R_full_stim_ph_trial, 3));
        assert(max(diffs(:) < 1e-4));

    elseif strcmp(phaseOEmode, 'oe')

        if odd(n_trials)
            R_full_stim_ph_trial(:,:,end-1) = mean(  R_full_stim_ph_trial(:,:,end-1:end), 3 );
            R_full_stim_ph_trial(:,:,end) = [];
            n_trials = n_trials -1 ;
        end
        idx_odd = 1:2:n_trials;
        idx_even = 2:2:n_trials;
        R_full_jackknifeTrials_stim_odd  = jackknifeAverageTrials(R_full_stim_ph_trial(:,:,idx_odd), 3);
        R_full_jackknifeTrials_stim_even = jackknifeAverageTrials(R_full_stim_ph_trial(:,:,idx_even), 3);


        R_full_jackknifeTrials = cell(n_stim, n_trials/2, 2);
        for jj = 1:n_stim
            R_full_jackknifeTrials(jj,:,1) = cellfun(@(x) x(jj,:), R_full_jackknifeTrials_stim_odd, 'un', 0);
            R_full_jackknifeTrials(jj,:,2) = cellfun(@(x) x(jj,:), R_full_jackknifeTrials_stim_even, 'un', 0);
        end

        R_av_jack = mean(cat(3, R_full_jackknifeTrials_stim_odd{:,:}, R_full_jackknifeTrials_stim_even{:,:}),3);
        diffs = abs(R_av_jack - mean(R_full_stim_ph_trial, 3));
        assert(max(diffs(:) < 1e-4));
    end
            
end
        

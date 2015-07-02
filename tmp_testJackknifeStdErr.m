% function tmp_testJackknifeStdErr


    n = 40;
    B = 10000;
    nTrials = 10;
    
    ccs = zeros(1,B);
    jack_stds = zeros(1,B);
    
    progressBar('init-', B, 30);
    for i = 1:B
        x1 = randn(nTrials, n);
        x2 = randn(nTrials, n);
        
        
        x1mean = mean(x1, 1);
        x2mean = mean(x2, 1);
        
        cc_i = pearsonR(x1mean, x2mean);
        
        
        x1_jacks = jackknifeAverageTrials( x1, 1 );
        x2_jacks = jackknifeAverageTrials( x2, 1 );

        x1_jacks2 = arrayfun(@(i) mean(  x1(setdiff(1:nTrials, i), : )  ,1), 1:nTrials, 'un', 0);
        x2_jacks2 = arrayfun(@(i) mean(  x2(setdiff(1:nTrials, i), : )  ,1), 1:nTrials, 'un', 0);
        
        assert(isequal(x1_jacks, x1_jacks2));
        assert(isequal(x2_jacks, x2_jacks2));
        
        ccs_jacks = zeros(1,nTrials);
        for j = 1:nTrials
            ccs_jacks(j) = pearsonR(x1_jacks{j}, x2_jacks{j});
        end
        jack_std = jackknifeStdErr(ccs_jacks, cc_i);
        
        ccs(i) = cc_i;
        jack_stds(i) = jack_std;

        
        
        progressBar(i);
    end
    
    3;











% end
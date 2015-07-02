function histChannelVariances

    Gids = getAllGids;
    nGids = length(Gids);
    
    chMeans = zeros(4,nGids);
    chCovs_after = zeros(4,nGids);    
    chCovs_before = zeros(4,nGids);    
    for i = 1:length(Gids)
        [chM, chCov] = getChannelMeansAndCovariance(Gids(i), 0);
        chCovs_after(:,i) = diag(chCov);        
        nm = getFileName('meansCov', Gids(i));
        nm = strrep(nm, 'MeanCov\', 'MeanCov\beforeFiltering\');
        S = load(nm);
        chCovs_before(:,i) = diag(S.channelCov);                
    end
    
    plotMeans = 0;
    plotStds = 1;
    
    if plotMeans
        figure(15);
        chMeans_abs_sm = sum(abs(chMeans),1);        
        hist(chMeans_abs_sm);
    end
    
    if plotStds
        chStds_before = sqrt(chCovs_before);
        chStds_after = sqrt(chCovs_after);
       hist((-chStds_before(:)+chStds_after(:))./(chStds_before(:))*100, 25); xlabel('% change')         
    end
    
%     figure(14);
%     hist(chCovs(1,:), 40);

    
%     hist(chMeans(:,:), 40);
    
%     hist(chMeans1(abs(chMeans1) < 10), 100);
    
    3
end
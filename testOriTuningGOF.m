function testOriTuningGOF
    global GC

    GC_good = [1648        1753        1753        2216        2344        2374        2374 2671        2771        2841        2903        3057        3085        4388  4470        4474        4488        4488        4488        4494        4494  4506        4506        4506        4522        4726        4796        4976  4998        5012        5048        5084        5158        5236        5248         5276
                   2           1           2           3           3           5           6 3           2           1           3           3           1           2  7           1           2           5           9           3           4 5           6           9           2           2           1           1 5           4           3           2           6           2           4  6]';
    GC_bad = [1735        2304        2344        2484              2661        2861     2883             4482        4718
           2           2           5           1                  1           2         4               1           4]';

    filename = getFileName('indiv', 'movie_fg');
    S_file = load(filename);        
%         varname = getName('celldata', Gids(i), cellIds(i));
%         v = S_file.(varname);        
%         tuningStats = calcDegreeOfTuningStats( v.OSP.R_full, v.OSP.bckgSamples, v.Gid);
        
    nGood = zeros(1,size(GC_good, 1));
    nBad  = zeros(1,size(GC_bad, 1));

%     s_good = cell(1,nGood);
    
    for i = 1:size(GC_good, 1)
        Gid = GC_good(i,1); cellId = GC_good(i,2);
        GC = [Gid, cellId];
        varname = getName('celldata', Gid, cellId);
        v = S_file.(varname);        
        tuningStats = calcDegreeOfTuningStats( v.OSP.R_full, v.OSP.bckgSamples, v.Gid);
                
        s_good(i) = tuningStats.oriStats_si;
        fprintf('*');
    end
       
%     s_bad = cell(1,nBad);
    for i = 1:size(GC_bad, 1);
        Gid = GC_bad(i,1); cellId = GC_bad(i,2);
        GC = [Gid, cellId];
        varname = getName('celldata', Gid, cellId);
        v = S_file.(varname);        
        tuningStats = calcDegreeOfTuningStats( v.OSP.R_full, v.OSP.bckgSamples, v.Gid);

        s_bad(i) = tuningStats.oriStats_si;
        fprintf('*');
    end
   
    otc_good = [s_good.otc_stats];
    otc_bad  = [s_bad.otc_stats];
    figure;
    
    gof_good = [otc_good.gof];
    gof_bad = [otc_bad.gof];
    h = hist2({realVals([otc_good.peakDist]), realVals([otc_bad.peakDist])}, [0:2:40]);
    figure;
    h = hist2({realVals([gof_good.adjrsquare]), realVals([gof_bad.adjrsquare])}, 20);
    h = hist2({realVals([gof_good.rmse]), realVals([gof_bad.rmse])}, 20);
    3;
    
       
       
       
       
       
end



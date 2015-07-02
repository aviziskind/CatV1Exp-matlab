for i = 1:length(Gids), 
    Gid = Gids(i);
%     i, 
%     try
%         detectChangesInInputCurrent(Gid); 
%     catch 
%         fprintf('Encountered an error on %d (Gid = %d)\n', i, Gid)
%         allErr_I = [allErr_I, i];
%         allErr_Gids = [allErr_Gids, Gid]; %#ok<*AGROW>
%     end

    sd = siteDataFor(Gid);
    dfInfo = sd.dataFileInfo;
    dfname_raw = dfInfo.dataFileName;
    pathname = rawDataDir(dfname_raw);
    
    dfname_raw = [pathname dfname_raw];
    matFileName = strrep(dfname_raw, '.dat', '.mat');

    %%        
    S_mat = load(matFileName);

    S = S_mat.currentChanges;
    P = S.plotData;
    
%     voltageChangeRange_nsamp = S.voltageChangeRange_nsamp;
%     shift_amt = voltageChangeRange_nsamp(1);
%     voltageChangeEnd_nsamp = diff(voltageChangeRange_nsamp);
%     
%     S_mat.currentChanges.shift_amt = shift_amt;
%     S_mat.currentChanges.voltageChangeEnd_nsamp = voltageChangeEnd_nsamp;
%     
%     save(matFileName, '-struct', 'S_mat');
%     
    allS(i) = S;

    3;

end


% allErr_Gids  [6276        6277        6301        6426        6476        6479]
% allErr_I =
%     62    63    76   104   110   113

function doAll
%     clear all;
%     global psthStatsSettings    
%     psthStatsSettings = struct('compressToList', false, 'compressToSparse', true, 'storeLogOfPvalues', true, ...
%         'storeUnassignedValuesAs', 'zeros', 'nudge0pval_amt', 1e-100, ...
%         'statsPrecision', 'double', 'ospPhCompressFcn', 'max');
% 
%     getPSTHwindowData;
%     
    
    
    clear all;
    global psthStatsSettings    
    psthStatsSettings = struct('compressToList', false, 'compressToSparse', true, 'storeLogOfPvalues', true, ...
        'storeUnassignedValuesAs', 'zeros', 'nudge0pval_amt', 1e-100, ...
        'statsPrecision', 'double', 'ospPhCompressFcn', 'mean');    

    getPSTHwindowData;

    


end
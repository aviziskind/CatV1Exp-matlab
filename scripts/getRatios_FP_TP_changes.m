function [allRatios, allFPchanges, allTPchanges, allfFPremain, allfTPremain] = getRatios_FP_TP_changes(resultStruct)
    allRatios = resultStruct.falsePos_fracRemoved ./ resultStruct.correctICspikes_fracRemoved;
    
%     allFPchanges = resultStruct.sorting_falsePosRate_change;
    a = resultStruct.fracOfClusterPruned;
    b = resultStruct.falsePos_fracRemoved;
    c = resultStruct.correctICspikes_fracRemoved;
    frfp = (b - a)./(1 - a);
%     assert(all ( frfp == allFPchanges ));
%     allFPchanges = frfp;
    allFPchanges = resultStruct.sorting_falsePosRate_change;


    frtp = (c - a)./(1 - a);
    allTPchanges = frtp;
    
    allfFPremain = 1-resultStruct.falsePos_fracRemoved;
    allfTPremain = 1-resultStruct.correctICspikes_fracRemoved;
%     medianRatio = median(allRatios);
    
%     isMedianRatioInf = isinf(medianRatio);
%     if isMedianRatioInf
%         medianRatio = max(allRatios(isfinite(allRatios)));
%         if isempty(medianRatio)
%             medianRatio = nan;
%         end
%     end
    %%
%     a = resultStruct.fracOfClusterPruned;
%     b = resultStruct.falsePos_fracRemoved;
%     c = resultStruct.correctICspikes_fracRemoved;
    
%     a = median(a);    b = median(b);    c = median(c);
    
%     fp_reduc2 = (b - a)./(1 - a);
    
%     allYields = resultStruct.falsePos_fracRemoved - resultStruct.fracOfClusterPruned;
%     allYields = resultStruct.sorting_falsePosRate_change ./ (resultStruct.fracOfClusterPruned .^2);
        
%     frfp = (b-a)./(1-a);
%     frtp = (c-a)./(1-a);
%     allYields = (1-a) ./ (1-c).^2;
    
%     medianYield = median(allYields);    
    
    
    
%     medianFPchange = mean(allFPchanges);
        
%     medianRatio = median(a);
%     medianFPchange = median(b);
%     medianYield = median(c);
    
end

%{

                 Gid_sort_prune: {1x24 cell}
                            Gid: [1x24 double]
           falsePos_fracRemoved: [1x24 double]
    correctICspikes_fracRemoved: [1x24 double]
            fracPruned_falsePos: [1x24 double]
            fracPruned_ICspikes: [1x24 double]
    sorting_falsePosRate_change: [1x24 double]
            fracOfClusterPruned: [1x24 double]
          fracOfClusterFalsePos: [1x24 double]
                nRefrSpikePairs: [7 7 15 3 58 6 4 8 4 4 3 5 3 2 7 5 18 3 7 22 23 5 4 5]
              pctRefrSpikePairs: [1x24 double]
              
  %}
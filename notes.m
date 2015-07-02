%%% Effectiveness of pruning:

% (1) Pruning:
printPruningResults
compareFalsePosNegBeforeAfterPruning
    
% (2) Cross pruning
printPruningResults;
compareFalsePosNegBeforeAfterCrossPruning;
    

%%% Error rates using different sorting features



%%% Filtering
compareFilteringMethods
testDifferentDetectionFilters

%%% Spike sorting: # of dimensions:
exploreBestNumberOfDimensions  % looking at actual dataa

measureEffectOfExtraDimensions  % simulations

%%% Spike sorting: all comparisons
compareSortingFeaturesPerformance
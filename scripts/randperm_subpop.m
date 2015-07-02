function all_idx = randperm_subpop(subpopulations, nUnits)
%     all_idx_orig = [subpopulations{:}];

    nPerSubpop = cellfun(@length, subpopulations);
    randperm_sub = arrayfun(@randperm, nPerSubpop, 'un', 0);
%     all_perms = cellfun(@(a,b) a(b), subpopulations, randperm_sub, 'un', 0);
%     nTotalCells = sum(nPerSubpop);
    highestIdx = max( cellfun(@max, subpopulations(nPerSubpop>0) ) );
    
    if ~exist('nUnits', 'var') || isempty(nUnits)
        nUnits = highestIdx;
    else
        assert(nUnits >= highestIdx);
    end
    
%     assert( all( cellfun(@(subpop) all(subpop < nUnits), subpopulations) ) );
    %%
    checkNoOverlap = 1;
    if checkNoOverlap
        used = false(1, nUnits);
        for i = 1:length(subpopulations);
            assert(~any(used(subpopulations{i})));
            used(subpopulations{i}) = true;
        end
    end
    
    all_idx = 1:nUnits;
    for i = 1:length(subpopulations);
        all_idx(subpopulations{i}) = all_idx(subpopulations{i}(randperm_sub{i}));
    end
    
    %check
    
%         for i = 1:length(subpopulations)
%             assert( isequal( sort(all_idx(subpopulations{i})), sort(subpopulations{i}) ));
%         end
    assert(isequal(sort(all_idx), 1:nUnits));    % this makes sure there are no overlaps, and no indices missing...
    
end

%{

[1 3 5], [2 4 6]

1 2 3 4 5 6
-->
5 4 6 1 2 3




%}

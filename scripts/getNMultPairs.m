function [combos_select, nApp] = getNMultPairs(nFet, nMult)
    combos = nchoosek(1:nFet, 2);
    ncombos = size(combos, 1);

    anyIsN = arrayfun(@(fi) any(combos == fi, 2), 1:nFet, 'un', 0); anyIsN = [anyIsN{:}];
    done = false;
    row_idx = 1;
    all_rows_idx = [];
%     curSum = zeros(ncombos, 1);
    while ~done

        sum_if_incl = sum(anyIsN([all_rows_idx, row_idx], :), 1);
        if all(sum_if_incl <= nMult)
            all_rows_idx = [all_rows_idx, row_idx]; %#ok<AGROW>
        end
        row_idx = row_idx+1;

        done = all(sum_if_incl == nMult) || (row_idx > ncombos);
    end

    combos_select = combos(all_rows_idx,:);
    nApp = sum(anyIsN(all_rows_idx, :), 1);
    
end
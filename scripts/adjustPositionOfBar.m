function y_bar_cent = adjustPositionOfBar(h, idx, newPos, barWidths, series_idx)
    3;
    h_patch = get(h, 'children');
    patch_v = get(h_patch, 'vertices');
    
    x = patch_v(:,1);
    y = patch_v(:,2);
    
%     idx_bars = cumsum( x(1:end-1) == 0 & (x(2:end) > 0));
%     i = find(idx_bars == idx,1);

    assert(length(idx) == length(newPos))
    patch_v_new = patch_v;
    y_new = y;
    nRowsPerBar = 5;
    
    cumBarWidths = [0, cumsum(barWidths)];
    barWidthTotal = sum(barWidths);
     
    for j = 1:length(idx)
        
        
        i = idx(j)*nRowsPerBar-3;

        %%
%         y_val = y(i);                       idx_y_val = find(y == y_val);
%         y_val_next = y(find(y > y_val, 1)); idx_y_val_next = find(y == y_val_next);
%         h_cur_cent = (y_val_next+y_val)/2;
%         y_shift = newPos(j) - h_cur_cent;
%     
%     
%         y_new(idx_y_val) = y_val + y_shift;
%         y_new(idx_y_val_next) = y_val_next + y_shift;
%         y_vals_rel = 
        
        y_val = y(i);                       idx_y_val = find(y == y_val);
        y_val_next = y(find(y > y_val, 1)); idx_y_val_next = find(y == y_val_next);
        
        h_new_top = newPos(j) - barWidthTotal/2 + cumBarWidths(series_idx); % + cumBarWidths(series_idx+1)/2;
        h_new_bot = newPos(j) - barWidthTotal/2 + cumBarWidths(series_idx+1);
        
%         h_cur_cent = (y_val_next+y_val)/2;
%         h_cur_top = (y_val_next+y_val)/2;
%         y_shift = newPos(j) - barWidthTotal/2 + cumBarWidths(series_idx) + cumBarWidths(series_idx+1)/2 - h_cur_cent;
        
        y_new(idx_y_val) = h_new_top;
        y_new(idx_y_val_next) = h_new_bot; %y_val_next + y_shift;

        y_bar_cent = (h_new_top + h_new_bot)/2;
    end
    patch_v_new(:,2) = y_new;
    set(h_patch, 'vertices', patch_v_new);
    
    
end
function drawLineThruBar(h)
    3;
    %%
    h_patch = get(h, 'children');
    patch_v = get(h_patch, 'vertices');
    
    x = patch_v(:,1);
    y = patch_v(:,2);
    
    idx_bars = cumsum( x(1:end-1) == 0 & (x(2:end) > 0));
    nbars = max(idx_bars);
    
    for bar_i = 1:nbars;            
        
        i = bar_i*5-3;

        
        y_val = y(i);                       idx_y_val = find(y == y_val);
        y_val_next = y(find(y > y_val, 1)); idx_y_val_next = find(y == y_val_next);    
        h_cur_cent = (y_val_next+y_val)/2;
        
        drawHorizontalLine(h_cur_cent);
    end
        
    
end
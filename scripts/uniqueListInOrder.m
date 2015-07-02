function [ux, idx] = uniqueListInOrder(x)
    [ux, idx] = uniqueList(x);
    [~, orig_idx] = sort( cellfun(@(is) is(1), idx) );
    ux = ux(orig_idx); 
    idx = idx(orig_idx);    
end

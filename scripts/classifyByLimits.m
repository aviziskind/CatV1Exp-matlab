function x_sections2 = classifyByLimits(x, limits)
    if ~issorted(limits)
        error('"Limits" input must be sorted');
    end
    if min((x(:)) < limits(1)) || max((x(:)) > limits(end))
        error('Values in x lie outside the categories');
    end
    x_orig = x;
    sizeX = size(x);
    x = x(:);
    nonnanInds = ~isnan(x);
    x = x(nonnanInds);
    
    small = 1e-3;
    limits(1)   = limits(1)  -small;
    limits(end) = limits(end)+small;
    
%     for i = 1:length(pvalSets)-1
%         inds = find( (pvalSets(i) <= cc_pval) & (cc_pval < pvalSets(i+1)) );
    x_sections = zeros(size(x));

    [x_sorted, inds_x] = sort(x);
        
    idx = elementsInRange(x_sorted, limits, 'index', '[)');
    for i = 1:length(idx)
        x_sections(  inds_x( idx{i} )  ) = i;
    end
    
    x_sections2 = x_orig;
    x_sections2(nonnanInds) = x_sections;
%     x_sections = resize(x_sections, sizeX);
    
    for i = 1:length(x)
        assert( ibetween( x(i), limits(x_sections(i)), limits(x_sections(i)+1) ));
    end


end



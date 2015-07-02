function tmp_recompressFiles(filename)

    S = load(filename);
    
    fn = fieldnames(S);
    s = whos('S');
    nbytes_before = s.bytes;

    progressBar('init-', length(fn))
    for i = 1:length(fn)
        v = S.(fn{i});        
        s2 = whos('v');
        n_before = s2.bytes;
        histVals1 = v.histVals;
        histVals2 = compress(single(decompress(histVals1)));
%         v.bckgSamples = compress(decompress(v.bckgSamples));
        v.histVals = histVals2;
        s2 = whos('v');
        n_after = s2.bytes;
        assert(n_after <= n_before+10);
%         if n_after <= n_before        
        S.(fn{i}) = v;
%         end
        progressBar;
    end
    
    s = whos('S');
    nbytes_after = s.bytes;
    3;
    save(filename, '-struct', 'S', '-v6');
    
end

%{
    
    

%}
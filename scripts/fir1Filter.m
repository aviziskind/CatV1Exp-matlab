function b = fir1Filter(n, Wn, ftype)
    % currently don't have the Signal Processing toolbox installed. this is
    % the result of calling 
    %    [b,a] = butter(n, Wn, ftype)
    % for the specific parameters that I use for filtering the datafile.
    %  (n = 1, Wn = 300/10000, ftype = 'high')

    if isequal({n, Wn, ftype}, {1, 300/10000, 'high'})
        b = [0.9549649940888687, -0.9549649940888687 ];
        a = [1,                  -0.90992998817773751];        
    elseif isequal({n, Wn, ftype}, {1, 300/5000, 'high'})
        b = [0.91363597298623778, -0.91363597298623778];
        a = [1,                   -0.82727194597247555];
        
    elseif isequal({n, Wn, ftype}, {1, [100 250]/10000, 'bandpass'})
        b = [0.023023722046118161, 0,  -0.023023722046118161];
        a = [1,  -1.9515420131006538,   0.95395255590776362];
        
    elseif isequal({n, Wn, ftype}, {1, [100 250]/5000, 'bandpass'})        
        b = [0.045035005911131229, 0, -0.045035005911131229];
        a = [1,  -1.9005056347917435, 0.909929988177737510];
    else         
        error('don''t have these filter parameters available');
    end

% [b,a] = butter(1, 300/10000, 'high');

end
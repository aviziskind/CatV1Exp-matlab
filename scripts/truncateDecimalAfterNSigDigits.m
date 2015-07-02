function str = truncateDecimalAfterNSigDigits(str, n)
    % truncate after n significant digits after decimal
    N = length(str);
    idx_dec = strfind(str, '.');

    i = idx_dec+1;
    if str2double(str) < 1
        while (i <= N) && str(i) == '0'
            i = i+1;
        end
    end
    if i <= N
        idx_first_sig_digit = i;
    else  % is just == 0;
        idx_first_sig_digit = idx_dec+1;
    end
        
    idx_nDigits = idx_first_sig_digit+n-1;    
    if idx_nDigits < N    
        if str2double(str(idx_nDigits+1)) >= 5 % have to round up.
            str(idx_nDigits+1) = '0'; % make sure doesn't round up in recursive call.
            valToAdd = 10^(idx_dec-idx_nDigits);            
            str = sprintf('%f', str2double(str)+valToAdd);
        end
        str = str(1:idx_nDigits);        
    end
end
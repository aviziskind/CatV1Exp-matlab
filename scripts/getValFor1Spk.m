function v = getValFor1Spk(r)
    u_r = unique(r(:));
    u_r_no0 = setdiff(u_r, 0);
    diffs = diff(u_r);
    v1 = min(diffs);
    
    if isempty(v1)
        v1 = 0;
    end
    
    v2 = v1;
    i = 1;
    done = false;    
    while ~done 
        v2 = v1/i;
        n_mult = u_r_no0/v2;    
        done = all( abs(n_mult-round(n_mult)) < 1e-5) || (i > 1000);
        if i > 1000
            warning('Exceeded threshold');
        end
        i = i + 1;
    end    
    
    if v1 > 0 && v1 ~= v2;
        fprintf('smallest diff = %.2g. gcd = %.2g\n', v1, v2)
%         keyboard;
    end
    v = v2;
end

% function g = gcd_all(v)
% 
%     g = v(1);
%     for i = 2:length(v)
%         g = gcd(g, v(i));
%     end
% end


%  4.4 6.6 9.9 
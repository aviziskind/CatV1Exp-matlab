function y2 = fermiSmooth(y1, k_hp, b_hp, circularFlag)
    % k_hp is the high-pass frequency.  for 60 phases, 30 does no smoothing
    % b_hp determines the steepness of the frequency drop-off
%     if (w == 0)
%         y2 = y1;
%         return;
%     end
    %%

    circularSmooth = exist('circularFlag', 'var') && ~isempty(circularFlag);
    N = length(y1);
%     n = max(ceil(5*w), 2); % how many std deviations of gaussian to actually implement

    if ~circularSmooth
        M = N*3;
        y1 = [ones(N,1)*y1(1); y1(:); ones(N,1)*y1(N)];
        k_hp = k_hp*3;
    else
        M = N;
    end
%     issym = @(x) all(x == conj(x(mod(length(x)-[1:length(x)]+1,length(x))+1)));
    ks = -M/2:M/2-1;
    F = fermi(ks, k_hp, b_hp);
        
    Y1 = fft(y1);
    y2 = ifft(Y1(:) .* F(:), 'symmetric');

%     figure(6);
%     xx = 1:length(y1); y1 = y1(:); y2 = y2(:);
%     plot(xx, y1, 'b.-', xx, y2, 'ro-')

    if ~circularSmooth
        y2 = y2([N+1:2*N]);        
    end
3;
    
    
end


function F = fermi(ks, k_hp, b_hp)

F = 1./(1+exp(-(k_hp-abs(ks))/b_hp));
F = fftshift(F);

end



%     if circularSmooth        
%         idx_pre  = mod([N-n+1:N]-1, N)+1;
%         idx_post = mod([1:n]-1, N)+1;  
%         y1 = [y1(idx_pre); y1; y1(idx_post)];
%     else
%         y1 = [y1(1)*ones(n, 1); y1; y1(end)*ones(n, 1)];
%     end
    
%     N = round(length(g)/2);
%     N2 = length(g)-N;
% 
%     f_conv = conv(y1, g);    
%     idx = [n+N:length(f_conv)-N2-n];
%     y2 = f_conv(idx);


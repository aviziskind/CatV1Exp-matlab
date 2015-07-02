function y = gaussOri_withMean(A, mu, sigma, B, x)
%     y = A .* exp( -(x - dir_pref_deg).^2 ./(2*(sigma).^2)) + B;
    y = A .* exp( -(x-mu).^2 ./(2*(sigma).^2)) + B;
end

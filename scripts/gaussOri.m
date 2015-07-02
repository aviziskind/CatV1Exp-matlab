function y = gaussOri(A, sigma, B, x)
%     y = A .* exp( -(x - dir_pref_deg).^2 ./(2*(sigma).^2)) + B;
    y = A .* exp( -(x ).^2 ./(2*(sigma).^2)) + B;
end

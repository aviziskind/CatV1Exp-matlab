function findAsciiNumber
    
    figure(1);
    k = 0;
    while ~isequal(k, 27) % Esc
        waitforbuttonpress;
        k = double(get(1, 'currentCharacter'));
        disp([num2str(k)]);
    end
    
    
end
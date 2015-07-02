allDids = zeros(size(allGids));
for i = 1:length(allGids)
    disp(outOf(i, length(allGids)));
    try
        allDids(i) = dbLookup('Did', 'Gid', allGids(i));    
    catch
        disp(['no Did for Gid = ' num2str(allGids(i))  ]);
        beep;
        pause(1);
    end
end


allSids = zeros(size(allGids));
for i = 1:allGids
    disp(outOf(i, length(allGids)));
    try
        allSids(i) = dbLookup('Sid', 'Gid', allGids(i));    
    catch
        Did = dbLookup('Did', 'Gid', allGids(i));  
        st = getStimulusTypeForDid(Did);
        disp(['none for Gid = ' num2str(allGids(i)) ' (' st ')' ]);
        beep;
        pause(1);
    end
end



% results (Did): all Groups have Datafiles.

% results: the following 9 Groups do not have Sids:
% Gid = 1561 (Grating)
% Gid = 1563 (Grating)
% Gid = 1565 (Grating)
% Gid = 1567 (Grating)
% Gid = 1571 (Grating)
% Gid = 1573 (Grating)
% Gid = 1575 (Grating)
% Gid = 1579 (Grating)
% Gid = 1581 (Grating)

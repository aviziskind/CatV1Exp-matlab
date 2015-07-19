function [timeWindow_out, timeWindow_str_out] = curTimeWindow(tw)
    
    persistent timeWindow

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curTimeWindowFile.mat'];
    numericTimeWindows = {[29, 62], [58, 91]};
    charTimeWindows = {'best', 'stimw', 'bestM2', 'bestM1', 'bestP1', 'bestP2'};
%     timeWindows = [numericTimeWindows, charTimeWindows];

    setType = (nargin == 1) && ~isempty(tw);
    
    if setType   % set current Comparison type 

        if ischar(tw)
            if ~any(strcmp(tw, charTimeWindows))
                error('Unknown time window')
            else
                timeWindow = tw;
            end
        elseif isnumeric(tw)
            idx = find(cellfun(@(t) isequal(t, tw), numericTimeWindows), 1);            
            if isempty(idx)
                if isequal(tw, [30 60])
                    tw = [29, 62];
                elseif isequal(tw, [60, 90])
                    tw = [58, 91];
                else                    
                    error('Invalid window');
                end
            end
            timeWindow = tw;
        else
            error('invalid timeWindow type:');
        end
        save(filename, 'timeWindow');
        
    else   % retrieve current Comparison type
        if isempty(timeWindow)
            if exist(filename, 'file')
                load(filename);
            else
                error('Time Window file does not exist');
            end
        end
        timeWindow_out = timeWindow;
        timeWindow_str_out = getTimeWindowStr(timeWindow); 
                
        if (nargin==1) && isempty(tw)
            timeWindow_out = timeWindow_str_out;
        end        

        
    end

end

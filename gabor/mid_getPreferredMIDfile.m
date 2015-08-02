function [mid_fileName, downSampFactor, frameMode] = mid_getPreferredMIDfile(Gid, cellId, timeWindow, trialMode, responseType)
    if nargin < 3
        timeWindow = 'best';
    end
    if nargin < 4
        trialMode = 'all'; % vs odd/even
    end
    if nargin < 5
        responseType = 'raw'; % vs odd/even
    end
    
    order_preference = {1  '1rep';
                        1, '2rep';                        
                        2  '1rep';
                        2  '2rep'};
        
    i=1; mid_fileName = '';
    while (i <= size(order_preference, 1)) && ~exist(mid_fileName, 'file');        
        downSampFactor = order_preference{i,1};
        frameMode      = order_preference{i,2};
        [mid_fileName, ~] = getName('MID_file', Gid, cellId, downSampFactor, frameMode, timeWindow, trialMode, responseType );
        i = i+1;        
    end
    if ~exist(mid_fileName, 'file')
        mid_fileName = '';
    end    
    
end
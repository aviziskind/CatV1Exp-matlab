function [p, z] = getValueProbFromDist(val, varargin)
%     p = getValueProbFromDist(val, nullDist)
%     p = getValueProbFromDist(val, nullDist_mean, nullDist_std)
%     p = getValueProbFromDist(..., tail)
%     [p, z] = getValueProbFromDist(...)
    
    if length(varargin{1}) > 1
        sampleDist = varargin{1};
        sample_mean = mean(sampleDist);
        sample_std  = std(sampleDist);
        varargin = varargin(2:end);
    else
        [sample_mean, sample_std] = varargin{1:2};
        varargin = varargin(3:end);
    end
    if ~isempty(varargin)
        tail = varargin{1};
    else
        tail = 'both';
    end

    erfc1 = @(x) .5*erfc(x/sqrt(2)); % integral of gaussian with variance 1 (instead of 1/2)        
            
    z = (val-sample_mean)/sample_std;
        
%     stat_pval_y = erfc1(val_y_norm);            
    
    switch tail
        case 'left',  p = erfc1(-z);
        case 'right', p = erfc1(z);
        case 'both',  p = 2*erfc1(abs(z)); 
        otherwise, error('unknown value for "tail"');
    end
    
    
end
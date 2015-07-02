function [spikeWaveforms, t_ms] = getSpikeWaveforms(Gid, cellId, normChannels, whitenChannels, matchDB, voltageFilterMode) %#ok<STOUT>
    
    % normChannels -- 0: do not normalize.  1: normalize channels separately.  2: normalize channels together 
    if ~exist('normChannels', 'var') || isempty(normChannels)
        normChannels = 0;
    elseif ischar(normChannels)
        normChannels = switchh(normChannels, {'norm', 'unnorm'}, [1, 0]);
    end    

    
    % whitenChannels -- 0: do not whiten channels.  1: whiten channels.  2: scale channels by std deviation
    if ~exist('whitenChannels', 'var') || isempty(whitenChannels) || isequal(whitenChannels, 0)
        whitenChannels  = 0;
        rescaleChannels = 0;
    elseif isnumeric(whitenChannels) && (whitenChannels == 1)
        whitenChannels = 1;
        rescaleChannels = 0;        
    elseif isnumeric(whitenChannels) && (whitenChannels == 2)        
        whitenChannels = 0;
        rescaleChannels = 1;
    elseif ischar(whitenChannels)        
        tfs = switchh(whitenChannels, {'raw', 'scl', 'ccw'}, {[0 0], [0 1], [1 0]});
        whitenChannels = tfs(1); rescaleChannels = tfs(2); 
    end
    
    if ~exist('matchDB', 'var') || isempty(matchDB)        
        matchDB = curMatchDB;
    end
    
    if ~exist('voltageFilterMode', 'var') || isempty(voltageFilterMode)
        voltageFilterMode = curVoltageFilter;
    end
    voltageFilterMode_str = getVoltageFilterName(voltageFilterMode, 'measure', Gid);
    
    file_opt = struct('voltageFilterMode_str', voltageFilterMode_str, 'matchDB', matchDB);
    
    spikeWaveformsFile  = getFileName('waveforms', Gid, 0, file_opt);
    load(spikeWaveformsFile);
    spikeWaveforms = single(spikeWaveforms); %#ok<NODEF>
    if exist('waveformFactor', 'var')
       spikeWaveforms = spikeWaveforms / waveformFactor;
    end
    [nT, nChannels, nSpk] = size(spikeWaveforms); 
    if exist('cellId', 'var') && ~isempty(cellId)
        [cellIds, cellIdxs] = getCellSorting(Gid, matchDB);
        id = find(cellId == cellIds,1);
        spikeWaveforms = spikeWaveforms(:,:, cellIdxs{id}); 
    end
           
    
    if whitenChannels 
        [channelMeans, channelCov, whitenMtx] = getChannelMeansAndCovariance(Gid, matchDB, voltageFilterMode);
                
        spikeWaveforms = reshape(permute(spikeWaveforms, [2 1 3]), [nChannels, nT*nSpk]);  % put channels in dim 1 so can matrix multiply.
        spikeWaveforms = whitenMtx * single(spikeWaveforms); % matrix conversion needs floating-point type
        
        spikeWaveforms = permute( reshape(spikeWaveforms, [nChannels, nT, nSpk]), [2 1 3]);  % back to original shape.                
        
    elseif rescaleChannels      
        [channelMeans, channelCov] = getChannelMeansAndCovariance(Gid, matchDB, voltageFilterMode);
        
        channelScaling = 1./sqrt(diag(channelCov));        
        spikeWaveforms = bsxfun(@times, channelScaling', spikeWaveforms);                       
    end
    
        
    if normChannels == 1   % normalize the spikes on each channel
        spikeWaveforms_norm = normV( spikeWaveforms, 1) / nT;
        spikeWaveforms_norm(spikeWaveforms_norm == 0) = 1;        
        spikeWaveforms = bsxfun(@rdivide, spikeWaveforms, spikeWaveforms_norm);                        


    elseif normChannels == 2  % normalize the entire (concatenated) spikes (usually after cross-channel whitening)
        spikeWaveforms = reshape(spikeWaveforms, [nT*nChannels, nSpk]);  

        spikeWaveforms_norm = normV( spikeWaveforms, 1) / (nT*nChannels);
        spikeWaveforms = bsxfun(@rdivide, spikeWaveforms, spikeWaveforms_norm);                                
        
        spikeWaveforms = reshape(spikeWaveforms, [nT, nChannels, nSpk]);          
    end
                
    
end


% 
% function [spikeWaveforms, t_ms] = getSpikeWaveForms(Gid, cellId, whitenChannels, matchDB) %#ok<STOUT>
%     if (nargin < 3) || isempty(matchDB)
%         matchDB = false;
%     end    
%     
%     if (nargin < 4) || isempty(whitenChannels) || isequal(logical(whitenChannels), 0);
%         wvfmType = 'raw';
%     elseif isequal(logical(whitenChannels), 1)
%         wvfmType = 'ccw';
%     elseif length(whitenChannels) > 1
%         wvfmType = 'both';
%     end        
%     
%     suffix = iff(matchDB, '_spkr', '');
%     subdir = iff(matchDB, 'spiker\', '');
%     spikeWaveformsFile  = [spikePath subdir 'Group_' num2str(Gid) '_spikeWaveforms' suffix '.mat'];
%     load(spikeWaveformsFile);
%     
%     if exist('cellId', 'var') && ~isempty(cellId)
%         [cellIds, cellIdxs] = getCellSorting(Gid, matchDB);
%         id = find(cellId == cellIds,1);
%         spikeWaveforms = spikeWaveforms(:,:, cellIdxs{id}); %#ok<NODEF>
%     end
%    
%     spikeWaveforms_raw = spikeWaveforms;
%     if any(strcmp(wvfmType, {'ccw', 'both'}))
%         % before we reshape the original spikeWaveforms, copy into the CCW matrix (and whiten the channels for computation below), 
%         [channelMeans, channelCov] = getChannelMeansAndCovariance(Gid);
%         [V,D] = eig(channelCov);
%         whitenMtx = sqrt(D)\V'; % whitening matrix        
%         
%         [nT, nChannels, nSpk] = size(spikeWaveforms);
%         spikeWaveforms_ccw = reshape(permute(spikeWaveforms, [2 1 3]), [nChannels, nT*nSpk]);  % put channels in dim 1 so can matrix multiply.
%         spikeWaveforms_ccw = whitenMtx * spikeWaveforms_ccw;
%         spikeWaveforms_ccw = permute( reshape(spikeWaveforms_ccw, [nChannels, nT, nSpk]), [2 1 3]);  % back to original shape.        
%     end
%     
%     switch wvfmType 
%         case 'raw',  spikeWaveforms = spikeWaveforms_raw;
%         case 'ccw',  spikeWaveforms = spikeWaveforms_ccw;
%         case 'both', spikeWaveforms = {spikeWaveforms_raw, spikeWaveforms_ccw};
%     end            
%     
%     
% end
% 

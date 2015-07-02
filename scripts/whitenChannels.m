function voltage_traces_whitened = whitenChannels(voltage_traces, whitenMtx)
                
    [nT, nChannels, nSpk] = size(voltage_traces);
    new_dim_order = [2 1 3];
    new_shape = [nChannels, nT*nSpk];

    back_to_orig_order = [2 1 3];
    orig_shape = [nChannels, nT, nSpk];
        % rearrange 
%         all_dims = 1:length(size(waveforms));
%         dims_reordered = [dim, setdiff(all_dims, dim)];
                
        
    voltage_traces_reshaped = reshape(permute(voltage_traces, new_dim_order), new_shape);

    if ~isfloat(voltage_traces_reshaped)
        voltage_traces_reshaped = single(voltage_traces_reshaped);
    end
    voltage_traces_whitened = whitenMtx * voltage_traces_reshaped; % matrix conversion needs floating-point type

    voltage_traces_whitened = permute( reshape(voltage_traces_whitened, orig_shape), back_to_orig_order); 

end


%         [channelMeans, channelCov, whitenMtx] = getChannelMeansAndCovariance(Gid);
%                 
%         spikeWaveforms = reshape(permute(spikeWaveforms, [2 1 3]), [nChannels, nT*nSpk]);  % put channels in dim 1 so can matrix multiply.
%         spikeWaveforms = whitenMtx * single(spikeWaveforms); % matrix conversion needs floating-point type
%         spikeWaveforms = permute( reshape(spikeWaveforms, [nChannels, nT, nSpk]), [2 1 3]);  % back to original shape.                

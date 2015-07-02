Gid = 4470;
[spikeWaveforms, t_ms] = getSpikeWaveforms(Gid);
spikeWaveforms1 = spikeWaveforms(:,:,randperm(15000));
[spikeWaveforms1a, t_ms_algn] = alignSpikeWaveforms(spikeWaveforms1, t_ms);

[coeffs1,  PCA_comps1, meanWaveform1] = calcWaveformCoefficients('pca', spikeWaveforms1, 5, 'separate');
[coeffs1a, PCA_comps1a, meanWaveform1a] = calcWaveformCoefficients('pca', spikeWaveforms1a, 5, 'separate');

% figure(1);
% subplot(4,8);
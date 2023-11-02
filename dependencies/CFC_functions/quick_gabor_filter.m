function ...
    filtered_LFP=quick_gabor_filter(raw_LFP,sampling_rate,center_frequency,fractional_bandwidth)

% function ...
%    filtered_LFP=quick_gabor_filter(raw_LFP,sampling_rate,center_frequency,fractional_bandwidth);
%
% input --
%     raw_LFP: 1xN real-valued time series of raw signal
%     sampling_rate: number of samples per second
%     center_frequency: mean frequency value of gaussian envelope of Gabor time-frequency atom
%     fractional_bandwidth: bandwidth parameter, frequency-domain standard deviation proportional to center_frequency*fractional_bandwidth
%
% output --
%     filtered_LFP: 1xN complex-valued time series of filtered signal

sp=get_signal_parameters('sampling_rate',sampling_rate,'number_points_time_domain',length(raw_LFP));

clear g
g.center_frequency=center_frequency; % Hz
g.fractional_bandwidth=fractional_bandwidth;
g.chirp_rate=0;
g=make_chirplet(...
    'chirplet_structure',g,...
    'signal_parameters',sp);

fs=filter_with_chirplet(...
    'raw_signal',raw_LFP,...
    'signal_parameters',sp,...
    'chirplet',g);

filtered_LFP=fs.time_domain;


% Load preprocessed data for each subject

% Then cut to pass the right time specification to ft_freqanalysis, which
% ignores toi specification in mtmfft
cfg = [];
cfg.toilim = [0, 1.5];
data_cut = ft_redefinetrial(cfg,data_tr_pp_ds);

% Compute fourier spectra for all pairs and trials
% fixed time window (t_ftimwin), single taper for all freqs
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foi        = 1:1:60; % changed from coh_pair_convol
cfg.taper      = 'hanning';
cfg.trials     = 'all';
cfg.pad        = 'nextpow2';
cfg.channel    = {'all'};
freq           = ft_freqanalysis(cfg, data_cut);

% Verify
cfg = [];
cfg.channel = 'ATS31';
cfg.parameter = 'powspctrm';
ft_singleplotER(cfg,freq)

% Connectivity analysis - all elecs at once since we're doing fft which is
% much lighter than mtmconvol
cfg            = [];
cfg.method     = 'coh'; % could be coh, ppc, or wpli_debiased
fdfourier      = ft_connectivityanalysis(cfg, freq);

% fdfourier.cohspctrm is chan_chan_freq - can be easily subset according to
% anatomical specification

fdfourier.epoch = epoch;
fdfourier.subset = subset;







cfg = [];
cfg.toilim = [0, 1.5];
data_cut = ft_redefinetrial(cfg,data_tr_pp_ds);
subset='ambiguous_gamb';

% time-frequency analysis
if strcmp(subset, 'all')
   trials = 'all'; % Passing a ones() vector causes problems with FT
   ntrials = sum(~bad_trials & ~timeout_trials);
   disp(['Analyzing ', num2str(ntrials), ' good trials']);
else
   trials = eval([subset]);
   trials = logical(reshape(trials, 1, []));
   %trials = trials(~bad_trials & ~timeout_trials);
   disp(['Analyzing ', num2str(sum(trials)), ' ', subset, ' trials out of ', num2str(sum(~bad_trials & ~timeout_trials)), ' good trials']);
end

cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foi        = 1:1:60; % changed from coh_pair_convol
cfg.taper      = 'hanning';
% cfg.t_ftimwin = 1.5; length of time window
% cfg.toi = 0.75; middle point
% cfg.t_ftimwin  = 0.25 .* ones(length(cfg.foi),1); % Width of analysis window - controls # lost data points at time edges
cfg.trials     = trials;
cfg.pad        = 'nextpow2';
cfg.channel    = 'all';
freq           = ft_freqanalysis(cfg, data_cut);

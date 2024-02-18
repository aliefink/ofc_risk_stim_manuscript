function [] = coh_pair_convol_LPFC_LJ(SUBID,elec_name,varargin)

%written by Lu Jin

%
% Calculate coherence between a specific electrode and all others in SUBID
% signal.mat: ft_definetrial -> ft_freqanalysis -> ft_connectivityanalysis
%
% Works on the ft_preprocessed dataset (ft_signal.mat), and is therefore
% subject to the same time windows
%
% Epoch is selected by using the appropriate ft_signal_epoch
% Event type is selected by subsetting trials - limited to binary classification: 'win', 'safebet', 'loss' (not xxx_ind)

disp(['Analyzing SUBID ', SUBID, ', elec ', elec_name]);

% Parse inputs
% Allows passing name-arg pairs, or using default values
defaultEpoch  = 'decision';  % rest, options, or decision
defaultSubset = 'all';      % all or binary (win_ind/loss_ind/safebet_ind etc.)
defaultMethod = 'fft';   % convol or fft
defaultTaper  = 'hanning';  % hanning or dpss
p = inputParser;

addRequired(p,'SUBID',@ischar);
addRequired(p,'elec_name',@ischar);
addParameter(p,'epoch', defaultEpoch, @ischar);
addParameter(p,'subset',defaultSubset,@ischar);
addParameter(p,'method',defaultMethod,@ischar);
addParameter(p,'taper',defaultTaper,@ischar);

parse(p,SUBID,elec_name,varargin{:});

epoch  = p.Results.epoch;
subset = p.Results.subset;
method = p.Results.method;
taper  = p.Results.taper;


% % Setup


% select data for epoch
%eval(['data = data_tr_',epoch]);
%eval(['data = signal_',epoch,';']);
%eval(['data = data_tr_pp_ds']);

% Load preprocessed data for each subject

% Then cut to pass the right time specification to ft_freqanalysis, which
% ignores toi specification in mtmfft

cfg = [];
toi = [-0.5,0]; % !!! set time window to cut the preprocessed data, 0.5 for options
%toi = [-1,0]; % !!! set time window to cut the preprocessed data, 1 for
%decision
cfg.toilim = toi; % set time window to cut the preprocessed data
data = ft_redefinetrial(cfg,data_tr_pp_ds);
eval(['data']);


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

% % Nasty patch to avoid a bad trial #136
% if strcmp(SUBID,'s32_IR55') && islogical(trials)
%     trials = trials(1:135);
% end


% Compute fourier spectra for all pairs and trials
% fixed time window (t_ftimwin), single taper for all freqs

if strcmp(taper,'hanning')
    disp('Using a hanning taper - optimal for f < 30Hz');
    % Computing fourier spectra
    % fixed time window (t_ftimwin), single taper for all freqs
    cfg            = [];
    cfg.output     = 'fourier';
    cfg.method     = ['mtm' method];
    cfg.foi        = 1:1:60; % changed from coh_pair_convol
    %cfg.foi        = create_freqs(1,30,1.5);
    %cfg.foi        = 10.^(log10(1):0.04:log10(30));
    cfg.taper      = 'hanning';
    % cfg.t_ftimwin = 1.5; length of time window
    % cfg.toi = 0.75; middle point
    % cfg.t_ftimwin  = 0.25 .* ones(length(cfg.foi),1); % Width of analysis window - controls # lost data points at time edges
    cfg.trials     = trials;
    cfg.pad        = 'nextpow2';
    cfg.channel    = {'all'};
    freq           = ft_freqanalysis(cfg, data);

elseif strcmp(taper,'dpss')
    disp('Using a DPSS multitaper - optimal for f > 30Hz');
    % Computing fourier spectra - multitapers
    % Not recommended for low fs - better for >30Hzx
    cfg            = [];
    cfg.output     = 'fourier';
    cfg.method     = ['mtm' method];
    cfg.taper      = 'dpss';
    cfg.trials     = trials;
    cfg.pad        = 'nextpow2';
    cfg.channel    = {'all'};
    cfg.foi        = 2:2:200; %used 2:2:200 for Decision epoch
%     cfg.toi        = toi; % need to be set
    t_ftimwin      = 5;
    cfg.t_ftimwin  = t_ftimwin./cfg.foi; % 5 cycles per time window
    cfg.tapsmofrq  = ((1/t_ftimwin) * cfg.foi)*2; % freq smoothing (x2)
    freq           = ft_freqanalysis(cfg, data);
        % disp(['Datalength=',num2str(cfg.t_ftimwin(i)),'s means min_smooth=',num2str(1/cfg.t_ftimwin(i)), '; foi=',num2str(cfg.foi(1)),'Hz, freq smoothing=',num2str(cfg.tapsmofrq(i))])
end
    
% Calculate coherence


if(~exist(out_dir))
    mkdir(out_dir);
end

% Connectivity analysis
cfg            = [];
cfg.method     = 'coh'; % could be coh, ppc, or wpli_debiased
cfg.channelcmb = {elec_name,data.label}; % Specify pairs of channels to analyze; wildcards may be used but both sets must be specified, e.g. {'HD_grid11*','HD_grid11*'}
fdfourier      = ft_connectivityanalysis(cfg, freq);
fdfourier.epoch = epoch;
fdfourier.subset = subset;
fdfourier.toi = toi;



% % Examine
% cfg  = [];
% cfg.parameter = 'cohspctrm';
% cfg.refchannel = 1;
% cfg.channel = 'all';
% cfg.title = elec_name;
% ft_singleplotTFR(cfg,fdfourier);

% Save
if strcmp(subset,'all')
    save([out_dir '/ft_',method,'_coh_',epoch,'_',elec_name{1},'.mat'],'fdfourier');
else
    save([out_dir '/ft_',method,'_coh_',epoch,'_',elec_name{1},'_',subset,'.mat'],'fdfourier');
end

end

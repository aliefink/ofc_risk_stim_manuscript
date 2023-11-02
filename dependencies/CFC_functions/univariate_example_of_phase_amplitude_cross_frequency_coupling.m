% example of phase-amplitude cross-frequency coupling

data_path='/Users/rcanolty/MATLAB/CFC_example/';
load GP1_gdat1_chan28_raw_signal

% useful signal parameters:
sp=get_signal_parameters('sampling_rate',2003,... % Hz
    'number_points_time_domain',length(raw_signal));

% time-frequency signal atom for low frequency filtering:
clear g1
g1.center_frequency=6; % Hz
g1.fractional_bandwidth=0.25;
g1.chirp_rate=0;
g1=make_chirplet(...
    'chirplet_structure',g1,...
    'signal_parameters',sp);

% time-frequency signal atom for high frequency filtering:
clear g2
g2.center_frequency=145; % Hz
g2.fractional_bandwidth=0.25;
g2.chirp_rate=0;
g2=make_chirplet(...
    'chirplet_structure',g2,...
    'signal_parameters',sp);

% filter raw signal at low frequency, extract phase:
fs=filter_with_chirplet(...
    'raw_signal',raw_signal,...
    'signal_parameters',sp,...
    'chirplet',g1);
lf_phase=angle(fs.time_domain);

% filter raw signal at high frequency, extract amplitude:
fs=filter_with_chirplet(...
    'raw_signal',raw_signal,...
    'signal_parameters',sp,...
    'chirplet',g2);
hfa=abs(fs.time_domain);

% filter high frequency amplitude time-series at low frequency, extract phase:
fs=filter_with_chirplet(...
    'raw_signal',hfa,...
    'signal_parameters',sp,...
    'chirplet',g1);
hf_phase=angle(fs.time_domain);

% make phase-triggered ERPs of raw signal and high-freq amplitude signals
% to see effect:

% find phases to trigger on:
[trash,maxima_indices]=find_maxima(lf_phase);

% epoch indices, +/- 1 second
einds=-sp.sampling_rate:sp.sampling_rate;

% drop phase indices that are outside bounds of averaging epoch:
maxima_indices(maxima_indices<=abs(einds(1)))=[];
maxima_indices(maxima_indices>=sp.number_points_time_domain-einds(end))=[];

% compute ERPs:
pt_erp=zeros(size(einds)); % (low freq) phase-triggered ERP of raw signal
pt_hfa_erp=zeros(size(einds)); % (low freq) phase-triggered ERP of high frequency amplitude
for e=1:length(maxima_indices)
    pt_erp=pt_erp+raw_signal(maxima_indices(e)+einds);
    pt_hfa_erp=pt_hfa_erp+hfa(maxima_indices(e)+einds);
end
pt_erp=pt_erp/length(maxima_indices);
pt_hfa_erp=pt_hfa_erp/length(maxima_indices);

figure;
h=subplot(1,2,1);
plot(einds,pt_erp);
xlim(einds([1 end]));
h=subplot(1,2,2);
plot(einds,pt_hfa_erp);
xlim(einds([1 end]));

% compute actual cross-frequency phase locking value (PLV):
actual_plv=abs(mean(exp(1i*(hf_phase-lf_phase))));

% compute an ensemble of surrogate PLVs to compare to actual value to
% establish statistical significance:
numsurrogate=10^3;
surrogate_plv=zeros(numsurrogate,1);
for s=1:numsurrogate
    disp(numsurrogate-s+1);
    shift=round(rand*sp.number_points_time_domain);
    surrogate_lf_phase=circshift(lf_phase,[0 shift]);
    surrogate_plv(s)=abs(mean(exp(1i*(hf_phase-surrogate_lf_phase))));
end
surrogate_plv=sort(surrogate_plv);

figure;
plot(surrogate_plv,'.k');
hold on;
plot(repmat(actual_plv,size(surrogate_plv)),'r');
hold off;


clear
clc

data_path='/Users/rcanolty/MATLAB/CFC_example/';

% for low frequency phases, use GP1, gdat1 channels: [7 8 15 16 24 35 36 49 62]
% lf_raw_signals=gdat([7 8 15 16 24 35 36 49 62],:);
%
% for high frequency channel, use GP1, gdat channel 16
% hf_raw_signal=gdat(16,:);

eval(['load ' data_path 'lf_raw_signals']);
eval(['load ' data_path 'hf_raw_signal']);

% useful signal parameters:
sp=get_signal_parameters('sampling_rate',2003,... % Hz
    'number_points_time_domain',length(hf_raw_signal));

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

% filter raw signal at high frequency, extract amplitude:
fs=filter_with_chirplet(...
    'raw_signal',hf_raw_signal,...
    'signal_parameters',sp,...
    'chirplet',g2);
hfa=abs(fs.time_domain);

% filter high frequency amplitude time-series at low frequency, extract phase:
fs=filter_with_chirplet(...
    'raw_signal',hfa,...
    'signal_parameters',sp,...
    'chirplet',g1);
hf_phase=angle(fs.time_domain);

% filter all lf_raw_signals at low frequency, extract phase:
lf_phase=zeros(size(lf_raw_signals));
for chan=1:size(lf_raw_signals,1)
    fs=filter_with_chirplet(...
        'raw_signal',lf_raw_signals(chan,:),...
        'signal_parameters',sp,...
        'chirplet',g1);
    lf_phase(chan,:)=angle(fs.time_domain);
end

all_phases=[hf_phase;lf_phase]; % note that high freq phases are first row

% fit multivariate model:
% second arguement is for sparsification; see Cadieu and Koepsell, 2010
% paper in Neural Computation
K = fit_model(all_phases,0.2*size(all_phases,1)); 
% K = a matrix of coupling parameters between all pairs of channels.  
% Kmn = kmn *exp[iumn]. kmn = coupling strength between channels m & n; 
% umn = preferred phase offset between channels m and n. 
%

mvals = K(1,2:end); % multivariate coupling of lf chans to hf chan

% fit several bivariate models:
bvals=zeros(size(mvals));
for chan=1:size(lf_raw_signals,1)
    bvals(chan)=mean(exp(1i*(hf_phase-lf_phase(chan,:))));%calculate bivariate PLVs.
end

figure;
for chan=1:size(lf_raw_signals,1)
    h=subplot(3,3,chan);
    [support,btemp]=make_von_mises_pdf(...
        angle(bvals(chan)),... von mises parameter mu
        abs(bvals(chan)),... von mises parameter kappa
        100);
    [support,mtemp]=make_von_mises_pdf(...
        angle(mvals(chan)),... von mises parameter mu
        abs(mvals(chan)),... von mises parameter kappa
        100);
    plot(support,btemp,'k'); % bivariate PLV-based distribution
    hold on;
    plot(support,mtemp,'r'); % multivariate PCE-based distribution
    hold off;
    xlim(pi*[-1 1]);
end



% fit surrogate versions of the multivariate model:
% this can take a while
% numsurrogate=200; % or more, more is better if you have time to run it
numsurrogate=20;
sur_K=zeros(size(K,1),size(K,2),numsurrogate);
sur_hf_phase=hf_phase;
ep=0.2*size(all_phases,1);
for s=1:numsurrogate
    disp(numsurrogate-s+1);
    sur_hf_phase=circshift(sur_hf_phase,[0 floor(rand*length(sur_hf_phase))]);
    sur_K(:,:,s)=fit_model([sur_hf_phase;lf_phase],ep);
end

pvalues=zeros(size(lf_raw_signals,1),1);
for chan=1:size(lf_raw_signals,1)
coupling_strength=abs(mvals(chan));
    sur_coupling_strength=flipud(sort(squeeze(abs(sur_K(1,1+chan,:)))));
    ind=find(coupling_strength<sur_coupling_strength,1,'last');
    if isempty(ind)
        pvalues(chan)=1/numsurrogate;
    else
        pvalues(chan)=ind/numsurrogate;
    end
end




    



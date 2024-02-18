function [sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count)


%Step One: Extract all normd plv matrices from all subjects
all_sig_plvs = zeros(40,19,sig_count); %amp x freq x elec
all_nonsig_plvs = zeros(40,19,nonsig_count); %amp x freq x elec

sig_num = 0; %tracker for num elecs - index for 3rd dimension of matrices above 
nonsig_num = 0; %tracker for num elecs - index for 3rd dimension of matrices above 

for s=1:length(subjs)
    %load cfc data for each subject
    subj_id = char(subjs{s});
    data = cfc_all_subj.(subj_id);
    %number of total electrodes
    num_sig_elecs = data.sig_elec_num; %num sig elecs for subj
    sig_elec_list = data.sig_elecs; %sig electrode ids
    
    for e=1:num_sig_elecs
        sig_num = sig_num+1;
        elec_id = sig_elec_list(e); %electrode id
        elec_idx = find(data.ofc_elecs==elec_id); %index of sig elec e in all data
        sig_norm_plv = data.norm_plvs{elec_idx};%get norm plv mat for sig elec
        all_sig_plvs(:,:,sig_num) = sig_norm_plv;
    end
    
    
    num_nonsig_elecs = data.nonsig_elec_num; %num sig elecs for subj
    nonsig_elec_list = data.nonsig_elecs; %sig electrode ids
    
    for e=1:num_nonsig_elecs
        nonsig_num = nonsig_num+1;
        elec_id = nonsig_elec_list(e); %electrode id
        elec_idx = find(data.ofc_elecs==elec_id); %index of sig elec e in all data
        nonsig_norm_plv = data.norm_plvs{elec_idx};%get norm plv mat for sig elec
        all_nonsig_plvs(:,:,nonsig_num) = nonsig_norm_plv;
    end
    
end

sig_mean_plv = mean(all_sig_plvs,3);
nonsig_mean_plv = mean(all_nonsig_plvs,3);

%%%%%%SIG
figure('Name','Significant Elec PLVs')
contourf(sig_mean_plv);

amp_f_array = 5:5:200;
phase_f_array = 2:20;

lim = max(sig_mean_plv,[],'all');
% Set locations of ticks and their labels.
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end), 'CLim', [0 lim] );
title(['Phase Amplitude Coupling of Significant Electrodes Across Subjects']);
xlabel('frequency for phase (Hz)');
ylabel('frequency for amplitude (Hz)');
h = colorbar;
set(get(h, 'ylabel'), 'string', 'PLV Zscore');


%%%%NONSIG
figure('Name','Not Significant Elec PLVs')

contourf(nonsig_mean_plv);

amp_f_array = 5:5:200;
phase_f_array = 2:20;


% Set locations of ticks and their labels.
lim = max(sig_mean_plv,[],'all');
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
    'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end), 'CLim', [0 lim] );
title(['Phase Amplitude Coupling of Non-Significant Electrodes Across Subjects']);
xlabel('frequency for phase (Hz)');
ylabel('frequency for amplitude (Hz)');
h = colorbar;
set(get(h, 'ylabel'), 'string', 'PLV Zscore');





end


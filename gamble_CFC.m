%% A. Fink CFC Analysis Fall 2023

% path loading
path_to_fieldtrip = '/Users/alexandrafink/Documents/MATLAB/SupportPackages/R2020b/fieldtrip-master/';
path_to_cfc ='/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/scripts/cfc/';
data_path ='/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/data/cfc_data/';
save_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/figs/cfc/';
fig_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/figs/cfc/results_methods/';


addpath(genpath(fig_path))
addpath(genpath(path_to_fieldtrip))
addpath(genpath(path_to_cfc))
addpath(genpath(data_path))
addpath(genpath(save_path))


%% Load CFC Results, Normalize Electrode-Level Data, + Calculate Pixel-Level Statistics 
    

%subj_ids in cfc dataset 
subjs = {'s02','s04','s06','s07','s08','s09','s10','s12','s13','s14','s15','s16','s17','s18','s19'};

%aggregate & normalize data + calculate individual pixel pvals
% cfc_all_subj = get_cfc_results_all(subjs,data_path); 
% save([data_path 'cfc_results/cfc_all_subj.mat'], 'cfc_all_subj')  %last saved 11/01/23 7:39PM
load([data_path 'cfc_results/cfc_all_subj.mat']) 

%% Group-Level Statistics 

% determine significant electrodes:
[sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'prop'); % can select pval 'zstat' or 'prop'
% save([data_path 'cfc_results/cfc_all_subj_sig_info.mat'], 'cfc_all_subj')%last saved 11/01/23 7:39PM

%0.95 cut off with corrected pvals only finds 22 sig elecs, up to 38 (or 33) after 
%removing amp correction, up to 41 or 36 if only max cluster is taken

% calculate & plot mean plv matrices for significant and nonsignificant electrodes
[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);



%% Extract mean plv data by freq band

[all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
    get_plv_freq_block_method(cfc_all_subj, subjs);


%% Get mean plv mat contours
%check bwlabeln and bwconncomp for cluster extraction

%extract contours from mean zscore sig elecs plot
cm = contourf(sig_mean_plv);
close all
contourTable = getContourLineCoordinates(cm);

% get indices from most sig clusters
max_z_all = max(contourTable.Level); %max level mean plv all sig elecs
sig_thresh = 2; %setting threshold for 2 std above mean - level > 2 *may need to raise (p < 0.05)
max_level = contourTable.Level > sig_thresh; %corresponds to p<0.001
max_levelX = contourTable.X(max_level); %extract x index of max contours (phase freq)
max_levelY = contourTable.Y(max_level); %extract y index of max contours (amp freq)

% get freq info from clust indices (must round contour index to whole number - actual data indices can't be exact contour curve)
sig_amp_idx = unique(round(max_levelY)); %indexes of contour amp freqs (rounded)
sig_phase_idx = unique(round(max_levelX)); %indexes of contour phase freqs (rounded)
amp_f_array = 5:5:200; % values of amp freqs
phase_f_array = 2:20; % values of phase freqs
sig_amps = amp_f_array(sig_amp_idx'); %amp freqs corresponding to amp index
sig_phases = phase_f_array(sig_phase_idx'); %phase freqs corresponding to phase index


%% Get contour indices by freq band
%define frequency bands
delta_theta = [2 8];
beta = [13 30];
gamma = [30 60];
hga = [60 200];
broadgamma = [30 200];


%% HGA x delta/theta PLV bar plots

%%%% extract the frequencies of amp/phase that are in the most significant
%%%% contours

cm = contourf(sig_mean_plv);
close all
contourTable = getContourLineCoordinates(cm);

% get indices from most sig clusters
max_z_all = max(contourTable.Level); %max level mean plv all sig elecs
sig_thresh = 2; %setting threshold for 2 std above mean - level > 2 *may need to raise (p < 0.05)
max_level = contourTable.Level > sig_thresh; %corresponds to p<0.001
max_levelX = contourTable.X(max_level); %extract x index of max contours (phase freq)
max_levelY = contourTable.Y(max_level); %extract y index of max contours (amp freq)

% get freq info from clust indices (must round contour index to whole number - actual data indices can't be exact contour curve)
sig_amp_idx = unique(round(max_levelY)); %indexes of contour amp freqs (rounded)
sig_phase_idx = unique(round(max_levelX)); %indexes of contour phase freqs (rounded)
amp_f_array = 5:5:200; % values of amp freqs
phase_f_array = 2:20; % values of phase freqs
sig_amps = amp_f_array(sig_amp_idx'); %amp freqs corresponding to amp index
sig_phases = phase_f_array(sig_phase_idx'); %phase freqs corresponding to phase index


[all_elec_plv_info,sig_elec_plv_info] = get_plv_freq_info_clust_method(cfc_all_subj, subjs, sig_amp_idx,sig_phase_idx);

% %define frequency bands
% delta_theta = [2 8];
% beta = [13 30];
% gamma = [30 60];
% hga = [60 200];
% 
% %find freqs by significant phases
% delta_theta_phases = sig_phases(sig_phases>= delta_theta(1)&sig_phases<= delta_theta(2));
% 
% %find freqs by significant amps
% beta_amps = sig_amps(sig_amps>= beta(1)&sig_amps< beta(2));
% gamma_amps = sig_amps(sig_amps>= gamma(1)&sig_amps< gamma(2));
% hga_amps = sig_amps(sig_amps>= hga(1)&sig_amps<= hga(2));
% 
% 
% 
% 
% 
% beta_delta_theta_amp_idx = find(ismember(amp_f_array,beta_amps));
% beta_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_beta_delta_theta_mean_cfc = sig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
% 
% gamma_delta_theta_amp_idx = find(ismember(amp_f_array,gamma_amps));
% gamma_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_gamma_delta_theta_mean_cfc = sig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
% 
% hga_delta_theta_amp_idx = find(ismember(amp_f_array,hga_amps));
% hga_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_hga_delta_theta_mean_cfc = sig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
% 


%% 



% statistical test between beta amp x delta-theta phase coupling
beta_delta_theta_amp_idx = find(ismember(amp_f_array,beta_amps));
beta_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
sig_beta_delta_theta_mean_cfc = sig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
nonsig_beta_delta_theta_mean_cfc = nonsig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);

[beta_delta_theta_clusters, beta_delta_theta_p_values, beta_delta_theta_t_sums, beta_delta_theta_permutation_distribution ] = permutest(sig_beta_delta_theta_mean_cfc,nonsig_beta_delta_theta_mean_cfc);
disp(beta_delta_theta_p_values)

mean_beta_delta_theta_cfc = mean(sig_beta_delta_theta_mean_cfc,'all');


% statistical test between gamma amp x delta-theta phase coupling
gamma_delta_theta_amp_idx = find(ismember(amp_f_array,gamma_amps));
gamma_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
sig_gamma_delta_theta_mean_cfc = sig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
nonsig_gamma_delta_theta_mean_cfc = nonsig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);

[gamma_delta_theta_clusters, gamma_delta_theta_p_values, gamma_delta_theta_t_sums, gamma_delta_theta_permutation_distribution ] = permutest(sig_gamma_delta_theta_mean_cfc,nonsig_gamma_delta_theta_mean_cfc);
disp(gamma_delta_theta_p_values)
mean_gamma_delta_theta_cfc = mean(sig_gamma_delta_theta_mean_cfc,'all');



% statistical test between hga amp x delta-theta phase coupling
hga_delta_theta_amp_idx = find(ismember(amp_f_array,hga_amps));
hga_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
sig_hga_delta_theta_mean_cfc = sig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
nonsig_hga_delta_theta_mean_cfc = nonsig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);

[hga_delta_theta_clusters, hga_delta_theta_p_values, hga_delta_theta_t_sums, hga_delta_theta_permutation_distribution ] = permutest(sig_hga_delta_theta_mean_cfc,nonsig_hga_delta_theta_mean_cfc);
disp(hga_delta_theta_p_values)
mean_hga_delta_theta_cfc = mean(sig_hga_delta_theta_mean_cfc,'all');

%% Figures single elec plv 
fig_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/figs/cfc/results_methods/';
addpath(genpath(fig_path))

data = cfc_all_subj.s06;

elec_num = 11;
plv_mat = data.norm_real_plvs{1,find(data.ofc_elecs==elec_num)};
plv_single_elec = figure('Name','zscore PLV one subj one elec')
contourf(plv_mat);
amp_f_array = 5:5:200;
phase_f_array = 2:20;
% Set locations of ticks and their labels.
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',14);
title(['Mean Phase Amplitude Coupling'],'FontSize',16);
xlabel('frequency for phase (Hz)','FontSize',14);
ylabel('frequency for amplitude (Hz)','FontSize',14);
h = colorbar;
set(get(h, 'ylabel'), 'string', 'PLV Zscore','FontSize',14);
saveas(plv_single_elec,[fig_path, num2str(elec_num), '_single_elec_plv.pdf'])



% plot mean amp ranges 
% 
% [sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, sig_elecs_by_subj, nonsig_elecs_by_subj, subjs, total_count,nonsig_count);


sig_plv_group = figure('Name','mean zscore PLV significant elecs across subj')
contourf(sig_mean_plv);
amp_f_array = 5:5:200;
phase_f_array = 2:20;
% Set locations of ticks and their labels.
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',14);
title(['Mean Phase Amplitude Coupling'],'FontSize',16);
xlabel('frequency for phase (Hz)','FontSize',14);
ylabel('frequency for amplitude (Hz)','FontSize',14);
h = colorbar;
set(get(h, 'ylabel'), 'string', 'PLV Zscore','FontSize',14);
saveas(sig_plv_group,[fig_path 'group_mean_plv_sig_elecs.pdf'])




%% Single electrode PLVs!


%% extract contours from contourf plot results 

[cm, h] = contourf(sig_mean_plv);
% contourTable = getContourLineCoordinates(cm);
[contourTable, contourArray] = getContourLineCoordinates(cm);

% get X,Y coords of most significant contours (> 3 std)
% max_level = contourTable.Level; == h.LevelList(7) | contourTable.Level == h.LevelList(8);
max_level = contourTable.Level > 3;
max_levelX = contourTable.X(max_level); %phases of most significant cluster
max_levelY = contourTable.Y(max_level); %amplitudes of most significant cluster

amp_f_array = 5:5:200;
phase_f_array = 2:20;

sig_phase_idx = unique(round(max_levelX));
sig_amp_idx = unique(round(max_levelY));
sig_phases = phase_f_array(sig_phase_idx');
sig_amps = amp_f_array(sig_amp_idx');

% Getting Mean PLV by freq band

%define frequency bands
delta_theta = [2 8];
alpha = [9 12];
beta = [13 29];
gamma = [30 60];
hga = [61 200];

%find freqs by significant phases
delta_theta_phases = sig_phases(sig_phases>= delta_theta(1)&sig_phases<= delta_theta(2));
alpha_phases = sig_phases(sig_phases>= alpha(1)&sig_phases<= alpha(2));

%find freqs by significant amps
alpha_amps = sig_amps(sig_amps>= alpha(1)&sig_amps<= alpha(2));
beta_amps = sig_amps(sig_amps>= beta(1)&sig_amps<= beta(2));
gamma_amps = sig_amps(sig_amps>= gamma(1)&sig_amps<= gamma(2));
hga_amps = sig_amps(sig_amps>= hga(1)&sig_amps<= hga(2));


% statistical test between beta amp x delta-theta phase coupling
beta_delta_theta_amp_idx = find(ismember(amp_f_array,beta_amps));
beta_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
sig_beta_delta_theta_mean_cfc = sig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
nonsig_beta_delta_theta_mean_cfc = nonsig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);

[beta_delta_theta_clusters, beta_delta_theta_p_values, beta_delta_theta_t_sums, beta_delta_theta_permutation_distribution ] = permutest(sig_beta_delta_theta_mean_cfc,nonsig_beta_delta_theta_mean_cfc);
disp(beta_delta_theta_p_values)

mean_beta_delta_theta_cfc = mean(sig_beta_delta_theta_mean_cfc,'all');


% statistical test between gamma amp x delta-theta phase coupling
gamma_delta_theta_amp_idx = find(ismember(amp_f_array,gamma_amps));
gamma_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
sig_gamma_delta_theta_mean_cfc = sig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
nonsig_gamma_delta_theta_mean_cfc = nonsig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);

[gamma_delta_theta_clusters, gamma_delta_theta_p_values, gamma_delta_theta_t_sums, gamma_delta_theta_permutation_distribution ] = permutest(sig_gamma_delta_theta_mean_cfc,nonsig_gamma_delta_theta_mean_cfc);
disp(gamma_delta_theta_p_values)
mean_gamma_delta_theta_cfc = mean(sig_gamma_delta_theta_mean_cfc,'all');



% statistical test between hga amp x delta-theta phase coupling
hga_delta_theta_amp_idx = find(ismember(amp_f_array,hga_amps));
hga_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
sig_hga_delta_theta_mean_cfc = sig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
nonsig_hga_delta_theta_mean_cfc = nonsig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);

[hga_delta_theta_clusters, hga_delta_theta_p_values, hga_delta_theta_t_sums, hga_delta_theta_permutation_distribution ] = permutest(sig_hga_delta_theta_mean_cfc,nonsig_hga_delta_theta_mean_cfc);
disp(hga_delta_theta_p_values)
mean_hga_delta_theta_cfc = mean(sig_hga_delta_theta_mean_cfc,'all');

%% Supplementary Figures 

% Histogram of permutation values 


surr = squeeze(cfc_all_subj.s02.perm_plvs{1,1}(1,1,:));

plv_zscore = cfc_all_subj.s02.norm_real_plvs{1,1}(1,1);

perm_fig = figure('Name','Single Electrode Single Pixel Permutation Distribution')
histogram(surr,20,'FaceColor','#D95319')
hold on 
xline(plv_zscore,'--k','Observed PLV zscore','LineWidth',1.5)
ylabel('Frequency')
xlabel('Surrogate Phase Locking Values')
title('Example Significant Permutation Distribution (n=1000)')
saveas(perm_fig,[fig_path 'permutation_example.pdf'])

% Another way to calculate pvalue - gives slightly different result...
% F = @(x)(exp (-0.5*(x.^2))./sqrt (2*pi));
% p = quad(F, plv_zscore, 100);





%% Plot average PLV of ALL electrodes within subj x risk pref 

risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
average_subj_plvs = [];

for s=1:length(subjs)
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    mean_plv_all_elec = mean(cell2mat(data.mean_norm_plvs));
    average_subj_plvs(1,s) = mean_plv_all_elec;
    

end 



lm = fitglm(average_subj_plvs,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(risk_prefs.',average_subj_plvs.','type','pearson');

scatter(average_subj_plvs,risk_prefs)
xlabel('Average PLV of All Electrodes by Subject')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))



%% Plot average PLV of Significant electrodes within subj x risk pref

risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
average_sig_plvs = [];


for s=1:length(subjs)
    subj_id = char(subjs(s));
    data = sig_elecs_by_subj{5,s}; %mean plv values for significant elecs
    mean_sig = mean(data);
    average_sig_plvs(1,s) = mean_sig;
    
end 


average_sig_plvs(isnan(average_sig_plvs)) = 0;

lm = fitglm(average_sig_plvs,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(risk_prefs.',average_sig_plvs.','type','pearson');

scatter(average_sig_plvs,risk_prefs)
xlabel('Average PLV of Significant Electrodes by Subject')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))



%% Plot average HGA-DELTA/THETA of Significant electrodes within subj x risk pref

risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
hga_cfc_bysubj = [];


for s=1:length(subjs)
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    hga_theta_cfc_means = [];
    for e=1:data.n_elecs
        plv = data.norm_real_plvs{1,e};
        hga_cfc_plv = plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
        hga_theta_cfc_means(1,e) = mean(hga_cfc_plv,'all');
    end 
    hga_cfc_bysubj(s) = mean(hga_theta_cfc_means);
end 



lm = fitglm(hga_cfc_bysubj,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(risk_prefs.',hga_cfc_bysubj.','type','pearson');

scatter(hga_cfc_bysubj,risk_prefs)
xlabel('Average HGA CFC PLV of Significant Electrodes by Subject')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))


%% Plot average GAMMA-DELTA/THETA of Significant electrodes within subj x risk pref

risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
gamma_cfc_bysubj = [];

for s=1:length(subjs)
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    gamma_theta_cfc_means = [];
    for e=1:data.n_elecs
        plv = data.norm_real_plvs{1,e};
        gamma_cfc_plv = plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
        gamma_theta_cfc_means(1,e) = mean(gamma_cfc_plv,'all');
    end 
    gamma_cfc_bysubj(s) = mean(gamma_theta_cfc_means);
end 



lm = fitglm(gamma_cfc_bysubj,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(risk_prefs.',gamma_cfc_bysubj.','type','pearson');

scatter(gamma_cfc_bysubj,risk_prefs)
xlabel('Average gamma CFC PLV of Significant Electrodes by Subject')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))



%% Plot average BETA-DELTA/THETA of Significant electrodes within subj x risk pref

risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
beta_cfc_bysubj = [];

for s=1:length(subjs)
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    beta_theta_cfc_means = [];
    for e=1:data.n_elecs
        plv = data.norm_real_plvs{1,e};
        beta_cfc_plv = plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
        beta_theta_cfc_means(1,e) = mean(beta_cfc_plv,'all');
    end 
    beta_cfc_bysubj(s) = mean(beta_theta_cfc_means);
end 



lm = fitglm(beta_cfc_bysubj,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(risk_prefs.',beta_cfc_bysubj.','type','pearson');

scatter(beta_cfc_bysubj,risk_prefs)
xlabel('Average beta CFC PLV of Significant Electrodes by Subject')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))




%% Plot average Significant cluster within subj x risk pref

risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
clust_cfc_bysubj = [];

for s=1:length(subjs)
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    clust_cfc_means = [];
    for e=1:data.n_elecs
        plv = data.norm_real_plvs{1,e};
        clust_cfc_plv = plv(sig_amp_idx',sig_phase_idx');
        clust_cfc_means(1,e) = mean(clust_cfc_plv,'all');
    end 
    clust_cfc_bysubj(s) = mean(clust_cfc_means);
end 



lm = fitglm(clust_cfc_bysubj,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(risk_prefs.',clust_cfc_bysubj.','type','pearson');

scatter(clust_cfc_bysubj,risk_prefs)
xlabel('Average Cluster CFC PLV of Significant Electrodes by Subject')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))



%% Average PLVs for NONsignificant electrodes

risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
average_nonsig_plvs = [];


for s=1:length(subjs)
    subj_id = char(subjs(s));
    data = nonsig_elecs_by_subj{4,s}; %mean plv values for significant elecs
    mean_nonsig = mean(data);
    average_nonsig_plvs(1,s) = mean_nonsig;
    
end 


average_nonsig_plvs(isnan(average_nonsig_plvs)) = 0;

lm = fitlm(average_nonsig_plvs,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(risk_prefs.',average_nonsig_plvs.','type','pearson');

scatter(average_nonsig_plvs,risk_prefs)
xlabel('Average PLV of Non Significant Electrodes by Subject')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))


%% Plot percent sig elecs within subj x risk pref 
risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% percent_sig_elecs = [0 34.69387755 85.71428571 9.090909091 70 42.85714286 33.33333333 0 33.33333333 50 33.33333333 0 0 0 0]; %prop
percent_sig_elecs = [0 30.6122449 85.71428571 9.090909091 70 14.28571429 0 0 33.33333333 50 33.33333333 0 0 0 0]; %tstat
% percent_sig_elecs = [20 44.89795918 57.14285714 72.72727273 50 42.857142860 0 33.33333333 50 33.33333333 0 0 16.66666667 0]; - gamma distributioncut offs


lm = fitlm(percent_sig_elecs,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(risk_prefs.',percent_sig_elecs.','type','pearson');

scatter(percent_sig_elecs,risk_prefs)
xlabel('Percent of Electrodes with Significant Coupling by Subject')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))







%% Step 2 find significant elecs:

[sig_elecs_by_subj,nonsig_elecs_by_subj,total_count,nonsig_count, cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'prop'); % p value calculations from proportion perms > pixel/num perms
%save([data_path 'cfc_results/cfc_all_subj.mat'], 'cfc_all_subj')  

% [sig_elecs_by_subj,nonsig_elecs_by_subj,total_count,nonsig_count] = find_sig_elecs(cfc_all_subj,subjs,'t_stat'); % p value calculations from t statistic not proprortion

%% Plot mean phase_amp_ranges 
    %by subject - all elecs, sig elecs, nonsig elecs
    %across subjects - all elecs, sig elecs, nonsig elecs
    
%inputs = normd plv matrix 
[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, sig_elecs_by_subj, nonsig_elecs_by_subj, subjs, total_count,nonsig_count);






%% Permutation Testing 

%permutation testing between mean plv matrices
% [mean_clusters, mean_p_values, mean_t_sums, mean_permutation_distribution ] = permutest(sig_mean_plv,nonsig_mean_plv);
% disp(mean_p_values);
% one_samp_plv_null = zeros(size(sig_mean_plv));
% [mean_clusters, mean_p_values, mean_t_sums, mean_permutation_distribution ] = permutest(sig_mean_plv,one_samp_plv_null);
% disp(mean_p_values);
% plv = cfc_all_subj.s02.real_plvs{1,1};
% perm_plv = mean(cfc_all_subj.s02.perm_plvs{1,1},3);
% [mean_clusters, mean_p_values, mean_t_sums, mean_permutation_distribution ] = permutest(plv,perm_plv);

%function [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(data,n_perm,tail,alpha_level,mu,reports,seed_state)

%permutation testing by specific frequency band


%% All Subjects CFC Calculations
%subjs = {'s02','s04','s06','s07','s08','s09','s10','s12','s13','s14','s15','s16','s17','s18','s19'};
%omit s13 briefly + s04 because num elecs
% subjs = {'s02','s06','s07','s08','s09','s10','s12','s13','s14','s15','s16','s17','s18','s19'};

subjs = {'s06','s07','s08','s09'};

num_perm = 1000;
cond_string = 'buttonpress_times';
start_window = -1000;
end_window = 0;

% parpool(3, 'IdleTimeout', Inf)

for s=1:length(subjs)
    subj_id = char(subjs(s));
    subj_id
    elec = load(strcat(data_path,subj_id,'/',subj_id,'_elec.mat'));  
    subj_globals = load(strcat(data_path,subj_id,'/',subj_id,'_subj_globals.mat'));

    anat_info = elec.elec.anat;
    ofc_idx = find(ismember(anat_info,'OrG')); %index of ofc ch 
    srate = subj_globals.srate; %this should now work for all subj
    
    parfor e=1:length(ofc_idx)
        elec_num = ofc_idx(e);
        [plv_matrix_real, plv_matrix_perm] = calculate_plv_real_perm_mats(data_path,...
            subj_id, cond_string, elec_num, srate, start_window, end_window, num_perm);
    end

end 

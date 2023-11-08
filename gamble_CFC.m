%% A. Fink CFC Analysis Fall 2023

% path loading
path_to_fieldtrip = '/Users/alexandrafink/Documents/MATLAB/SupportPackages/R2020b/fieldtrip-master/';
path_to_scripts ='/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/scripts/cfc/';
data_path ='/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/data/cfc_data/';
save_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/figs/cfc/';
fig_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/figs/cfc/results_methods/';



addpath(genpath(path_to_fieldtrip))
addpath(genpath(path_to_scripts))
addpath(genpath(data_path))
addpath(genpath(save_path))
addpath(genpath(fig_path))


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


% extract significant pixels:
% [sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_pixels(cfc_all_subj,subjs,'prop'); % can select pval 'zstat' or 'prop'
% [sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);



%% Extract mean plv data by freq band

[all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
    get_plv_freq_block_method(cfc_all_subj, subjs);













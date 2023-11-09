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
% save([data_path 'cfc_results/cfc_all_subj.mat'], 'cfc_all_subj')  %last saved 11/08/23 7:04PM
load([data_path 'cfc_results/cfc_all_subj.mat']) 

%% Group-Level Statistics 

%%%%%CLUSTER METHOD
% determine significant electrodes:
    %pval_type can be either 'prop' or 'zstat' - for paper analysis prop
    %plot_elec can either be 'yes' or 'no' - if yes will plot norm plv contourf & save
[sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'prop','no'); % can select pval 'zstat' or 'prop'
%run with saving elec figs - (last ran 11/8/23 7:09PM)
% [sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'prop','yes'); % can select pval 'zstat' or 'prop'
% save([data_path 'cfc_results/cfc_all_subj_sig_info.mat'], 'cfc_all_subj')%last saved 11/08/23 8:29PM
% save([data_path 'cfc_results/sig_elecs_by_subj.mat'], 'sig_elecs_by_subj')%last saved 11/08/23 8:29PM
% save([data_path 'cfc_results/nonsig_elecs_by_subj.mat'], 'nonsig_elecs_by_subj')%last saved 11/08/23 8:29PM
% save([data_path 'cfc_results/find_sig_elecs_header.mat'], 'find_sig_elecs_header')%last saved 11/08/23 8:29PM

%%%%%PIXEL METHOD
%extract significant pixels: - gives 121/130 electrodes as significant (11/8/23 8:15PM)
% [sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_pixels(cfc_all_subj,subjs,'prop'); % can select pval 'zstat' or 'prop'
% [sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);

% calculate & plot mean plv matrices for significant and nonsignificant electrodes
[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);


%% Extract mean plv data by freq band


[all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
    get_plv_freq_block_method(cfc_all_subj, subjs);

% save([data_path 'cfc_results/all_elec_plv_info.mat'], 'all_elec_plv_info')
% save([data_path 'cfc_results/sig_elec_plv_info.mat'], 'sig_elec_plv_info')
% save([data_path 'cfc_results/nonsig_elec_plv_info.mat'], 'nonsig_elec_plv_info')
% save([data_path 'cfc_results/phase_freqs.mat'], 'phase_freqs')
% save([data_path 'cfc_results/amp_freqs.mat'], 'amp_freqs')






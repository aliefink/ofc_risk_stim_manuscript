%% A. Fink CFC Analysis Fall 2023

%step-by-step analysis pipeline for risk project cross frequency coupling
%analysis 


%% Directory setup 


% path loading - update paths to local
path_to_fieldtrip = '/Users/alexandrafink/Documents/MATLAB/SupportPackages/R2020b/fieldtrip-master/'; %MATLAB fieldtrip package
data_path ='/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/data/cfc_data/'; %raw signal data
save_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/figs/cfc/'; %fig saving path
fig_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/figs/cfc/results_methods/'; %fig saving path


addpath(genpath(fig_path))
addpath(genpath(path_to_fieldtrip))
addpath(genpath(data_path))
addpath(genpath(save_path))


%% Step One: Calculate CFC + Surrogate CFC Matrices
%longest step in pipeline (# permutations is rate limiting step)
%main dependencies - CFC_functions package, calculate_plv_real_perm_mats

subjs = {'s02','s04','s06','s07','s08','s09','s10','s12','s13','s14','s15','s16','s17','s18','s19'};


num_perm = 1000; %200 will be faster
cond_string = 'buttonpress_times';
start_window = -1000;
end_window = 0;

%parpool(3, 'IdleTimeout', Inf) %paralellization option to speed up computation - parallelize by elec within subj

for s=1:length(subjs)
    subj_id = char(subjs(s));
    subj_id %print subj_id to track loop progress
    elec = load(strcat(data_path,subj_id,'/',subj_id,'_elec.mat')); %electrode info for each subj
    subj_globals = load(strcat(data_path,subj_id,'/',subj_id,'_subj_globals.mat')); %recording info for each subj

    anat_info = elec.elec.anat; %anatomical localization info for each elec
    ofc_idx = find(ismember(anat_info,'OrG')); %index of ofc ch 
    srate = subj_globals.srate; %sampling rate of raw signal data
    
    parfor e=1:length(ofc_idx)
        elec_num = ofc_idx(e); %electrode id number from elec.mat
        elec_num %print elec_num to track loop progress 
        [plv_matrix_real, plv_matrix_perm] = calculate_plv_real_perm_mats(data_path,...
            subj_id, cond_string, elec_num, srate, start_window, end_window, num_perm); %PLV + surrogate PLV calculations
            %function autosaves real + permutation PLV matrices for each electrode
    end

end 

%% Step Two: Calculate CFC Results + Aggregrate Across Subjects
%for each subject: normalized PLV matrix, pixel-level p value
%calculations, mean PLVs within and across electrodes

%main dependencies: fdr_bh package, get_cfc_results_all

cfc_all_subj = get_cfc_results_all(subjs,data_path);
% save([data_path 'cfc_results/cfc_all_subj.mat'], 'cfc_all_subj') 


%% Step Three: Identify, Visualize, & Summarize Electrodes with Significant CFC
%significance can be determined from two options of pixel-level p value calculations -
%prop = (num_perm values > real pixel)/total perm ; 
%zstat = p value from z statistic (normalized pixel)

%main dependencies: fdr_bh package, getContourLineCoordinates package, 
    %find_sig_elecs, plot_mean_phase_amp_ranges, get_plv_freq_block_method

[sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, ...
    find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'zstat'); % can select pval 'zstat' or 'prop'

% calculate & plot mean plv matrices for significant and nonsignificant electrodes
[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = ...
    plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);

%extract mean PLV data by frequency band for all electrodes + sig/nonsig
[all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
    get_plv_freq_block_method(cfc_all_subj, subjs); %can use block or cluster method


%% Step Four: 







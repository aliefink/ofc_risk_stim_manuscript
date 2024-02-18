function [all_elec_plv_info,sig_elec_plv_info] = get_delta_theta_gamma_clust_info(cfc_all_subj, subjs, sig_amp_idx,sig_phase_idx)
%outputs: all_elec_plv_info, sig_elec_plv_info

% frequency-specific information - define for all subj
amp_f_array = 5:5:200; % values of amp freqs
phase_f_array = 2:20; % values of phase freqs
%define frequency bands of interest 
delta_theta = [2 8]; % phase freq of interest
beta = [13 30];
gamma = [30 70];
hga = [70 200];
broadgamma = [30 200];
%define phase freqs of interest (delta_theta only)
delta_theta_phase_idx = sig_phase_idx(phase_f_array(sig_phase_idx')<=delta_theta(2)); %all delta_theta indices of interest
%define amp freqs of interest (delta_theta only)
beta_amp_idx = sig_amp_idx(amp_f_array(sig_amp_idx')>=beta(1) & amp_f_array(sig_amp_idx')<beta(2)); %beta amp idx 
gamma_amp_idx = sig_amp_idx(amp_f_array(sig_amp_idx')>=gamma(1) & amp_f_array(sig_amp_idx')<gamma(2)); %gamma amp idx
hga_amp_idx = sig_amp_idx(amp_f_array(sig_amp_idx')>=hga(1)& amp_f_array(sig_amp_idx')<=hga(2)); %hga 
broadgamma_amp_idx = sig_amp_idx(amp_f_array(sig_amp_idx')>=gamma(1)& amp_f_array(sig_amp_idx')<=hga(2)); %gamma+hga combo


% means from entire comodulograms
all_elec_mean_plv_raw = []; %mean raw PLV for all electrodes within subj
all_elec_mean_plv_zscore = []; %mean zscore PLV for all electrodes within subj
sig_elec_mean_plv_zscore = []; %mean zscore PLV for significant electrodes within subj
nonsig_elec_mean_plv_zscore = []; %mean zscore PLV for NOT significant electrodes within subj

% means from specific frequency bands
all_elec_mean_total_cluster = [];
all_elec_mean_deltatheta_beta = [];
all_elec_mean_deltatheta_gamma = [];
all_elec_mean_deltatheta_hga = [];
all_elec_mean_deltatheta_broadgamma = [];

sig_elec_mean_total_cluster = [];
sig_elec_mean_deltatheta_beta = [];
sig_elec_mean_deltatheta_gamma = [];
sig_elec_mean_deltatheta_hga = [];
sig_elec_mean_deltatheta_broadgamma  = [];
sig_subjs = {}; %array of subj with at least 1 significant electrode
sig_subj_count = 0;

for s=1:length(subjs) %iterate through all subjects   
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    
    %extract mean all elec info 
    all_elec_mean_plv_raw(s) = data.mean_real_plv_all_elec;
    all_elec_mean_plv_zscore(s) = data.mean_norm_plv_all_elec;
    sig_elec_mean_plv_zscore(s) = data.sig_mean_norm_plv; % mean sig elec norm plv
    nonsig_elec_mean_plv_zscore(s) = data.nonsig_mean_norm_plv; % mean nonsig elec norm plv
    
    subj_total_cluster_mean = [];
    subj_beta_cluster_mean = [];
    subj_gamma_cluster_mean = [];
    subj_hga_cluster_mean = [];
    subj_broadgamma_cluster_mean = [];
    sig_total_cluster_mean = [];
    sig_beta_cluster_mean = [];
    sig_gamma_cluster_mean = [];
    sig_hga_cluster_mean = [];
    sig_broadgamma_cluster_mean = [];
    
    for e=1:data.n_elecs
        elec = data.ofc_elecs(e);
        plv = data.norm_plvs{e}; %electrode plv is normalized plv matrix 
        total_cluster_mat = plv(sig_amp_idx,sig_phase_idx); %extract plv info from mean subj cluster inputs
        total_cluster_mean = mean(total_cluster_mat,'all'); %mean plv of total cluster 
        subj_total_cluster_mean(e) = total_cluster_mean; %add to subj mean vector
        
        %define clusters based on specific freq indices defined above 
        %beta
        beta_cluster_mat = plv(beta_amp_idx,delta_theta_phase_idx);
        beta_cluster_mean = mean(beta_cluster_mat,'all'); %mean plv of total cluster 
        subj_beta_cluster_mean(e) = beta_cluster_mean; %add to subj mean vector
        %gamma
        gamma_cluster_mat = plv(gamma_amp_idx,delta_theta_phase_idx);
        gamma_cluster_mean = mean(gamma_cluster_mat,'all'); %mean plv of total cluster 
        subj_gamma_cluster_mean(e) = gamma_cluster_mean; %add to subj mean vector
        %hga
        hga_cluster_mat = plv(hga_amp_idx,delta_theta_phase_idx);
        hga_cluster_mean = mean(hga_cluster_mat,'all'); %mean plv of total cluster 
        subj_hga_cluster_mean(e) = hga_cluster_mean; %add to subj mean vector
        %broadgamma
        broadgamma_cluster_mat = plv(broadgamma_amp_idx,delta_theta_phase_idx);
        broadgamma_cluster_mean = mean(broadgamma_cluster_mat,'all'); %mean plv of total cluster 
        subj_broadgamma_cluster_mean(e) = broadgamma_cluster_mean; %add to subj mean vector
        
        if sum(find(data.sig_elecs==elec)) ~= 0 %check for elec in sig elecs index - if it is there, sum =1, if not =0
            sig_idx = find(data.sig_elecs==elec);
            sig_beta_cluster_mean(sig_idx) = beta_cluster_mean;
            sig_gamma_cluster_mean(sig_idx) = gamma_cluster_mean; %add to subj mean vector
            sig_hga_cluster_mean(sig_idx) = hga_cluster_mean;
            sig_broadgamma_cluster_mean(sig_idx) = broadgamma_cluster_mean;
            sig_total_cluster_mean(e) = total_cluster_mean; %add to subj mean vector
        
        end 

    end %elec loop
        

    % means from specific frequency bands
    all_elec_mean_total_cluster(s) = mean(subj_total_cluster_mean);
    all_elec_mean_deltatheta_beta(s) = mean(subj_beta_cluster_mean);
    all_elec_mean_deltatheta_gamma(s) = mean(subj_gamma_cluster_mean);
    all_elec_mean_deltatheta_hga(s) = mean(subj_hga_cluster_mean);
    all_elec_mean_deltatheta_broadgamma(s) = mean(subj_broadgamma_cluster_mean);

    sig_elec_mean_total_cluster(s) = mean(sig_total_cluster_mean);
    sig_elec_mean_deltatheta_hga(s) = mean(sig_beta_cluster_mean);
    sig_elec_mean_deltatheta_gamma(s) = mean(sig_gamma_cluster_mean);
    sig_elec_mean_deltatheta_beta(s) = mean(sig_hga_cluster_mean);
    sig_elec_mean_deltatheta_broadgamma(s) = mean(sig_broadgamma_cluster_mean);
    if mean(sig_total_cluster_mean) ~= 0
         sig_subj_count = sig_subj_count+1;
         sig_subjs{sig_subj_count} = subj_id; 
    end 

        
end %subj loop

    %add cluster indexing info to all_elec_plv_info
    all_elec_plv_info.subjs = subjs;
    all_elec_plv_info.sig_amp_idx = sig_amp_idx;
    all_elec_plv_info.sig_phase_idx =  sig_phase_idx;
    all_elec_plv_info.beta_amp_idx = beta_amp_idx;
    all_elec_plv_info.gamma_amp_idx = gamma_amp_idx;
    all_elec_plv_info.hga_amp_idx = hga_amp_idx;
    all_elec_plv_info.broadgamma_amp_idx = broadgamma_amp_idx;
    %add all elec all subj info to all_elec_plv_info
    all_elec_plv_info.mean_plv_raw = all_elec_mean_plv_raw; %mean raw PLV for all electrodes within subj
    all_elec_plv_info.mean_plv_norm = all_elec_mean_plv_zscore; %mean zscore PLV for all electrodes within subj
    all_elec_plv_info.mean_total_cluster = all_elec_mean_total_cluster;
    all_elec_plv_info.mean_deltatheta_beta = all_elec_mean_deltatheta_beta;
    all_elec_plv_info.mean_deltatheta_gamma = all_elec_mean_deltatheta_gamma;
    all_elec_plv_info.mean_deltatheta_hga = all_elec_mean_deltatheta_hga;
    all_elec_plv_info.mean_deltatheta_broadgamma = all_elec_mean_deltatheta_broadgamma;
    
    %add cluster indexing info to sig_elec_plv_info
    sig_elec_plv_info.sig_subjs = sig_subjs;
    sig_elec_plv_info.sig_amp_idx = sig_amp_idx;
    sig_elec_plv_info.sig_phase_idx =  sig_phase_idx;
    sig_elec_plv_info.beta_amp_idx = beta_amp_idx;
    sig_elec_plv_info.gamma_amp_idx = gamma_amp_idx;
    sig_elec_plv_info.hga_amp_idx = hga_amp_idx;
    sig_elec_plv_info.broadgamma_amp_idx = broadgamma_amp_idx;
    %add sig elec info to sig_elec_plv_info
    sig_elec_plv_info.mean_plv_norm = sig_elec_mean_plv_zscore; %mean zscore PLV for all electrodes within subj
    sig_elec_plv_info.mean_total_cluster = sig_elec_mean_total_cluster;
    sig_elec_plv_info.mean_deltatheta_beta = all_elec_mean_deltatheta_beta;
    sig_elec_plv_info.mean_deltatheta_gamma = sig_elec_mean_deltatheta_gamma;
    sig_elec_plv_info.mean_deltatheta_hga = all_elec_mean_deltatheta_hga;
    sig_elec_plv_info.mean_deltatheta_broadgamma = sig_elec_mean_deltatheta_broadgamma;
    

end 
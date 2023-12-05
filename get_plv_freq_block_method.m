function [all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = get_plv_freq_block_method(cfc_all_subj, subjs)

%%%extract phase freq information for bands of interest (delta_theta alpha beta)

% phase frequencies of comodulogram 
phase_f_array = 2:20; % values of phase freqs
%define phase freqs of interest (delta_theta alpha beta)
delta_theta = [2 9]; % phase freq of interest
alpha = [9 13];
beta = [13 30];
%phase freq indices 
delta_theta_idx = find((phase_f_array>=delta_theta(1)&phase_f_array<delta_theta(2)));
alpha_idx = find((phase_f_array>=alpha(1)&phase_f_array<alpha(2)));
beta_idx = find((phase_f_array>=beta(1)&phase_f_array<beta(2)));
%group phase freqs into cell array
phase_freqs = {}; %first row is freq band name, second is freq band idx, third is freq band values
%delta-theta
phase_freqs(1,1) = {'delta_theta'};
phase_freqs(2,1) = {delta_theta_idx};
phase_freqs(3,1) = {phase_f_array((phase_f_array>=delta_theta(1)&phase_f_array<delta_theta(2)))};
%alpha
phase_freqs(1,2) = {'alpha'};
phase_freqs(2,2) = {alpha_idx};
phase_freqs(3,2) = {phase_f_array((phase_f_array>=alpha(1)&phase_f_array<alpha(2)))};
%beta
phase_freqs(1,3) = {'beta'};
phase_freqs(2,3) = {beta_idx};
phase_freqs(3,3) = {phase_f_array((phase_f_array>=beta(1)&phase_f_array<beta(2)))};

%%%extract amp freq information for bands of interest (gamma,hga,broadgamma)
% amp frequencies of comodulogram 
amp_f_array = 5:5:200; % values of amp freqs
%define amp freqs of interest (gamma,hga,broadgamma)
gamma = [30 70];
hga = [70 200];
broadgamma = [30 200];
%amp freq indices
gamma_idx = find((amp_f_array>=gamma(1)&amp_f_array<=gamma(2)));
hga_idx = find((amp_f_array>hga(1)&amp_f_array<=hga(2)));
broadgamma_idx = find((amp_f_array>=broadgamma(1)&amp_f_array<=broadgamma(2)));
%group amp freqs into cell array
amp_freqs = {}; %first row is freq band name, second is freq band idx, third is freq band values
%delta-theta
amp_freqs(1,1) = {'gamma'};
amp_freqs(2,1) = {gamma_idx};
amp_freqs(3,1) = {amp_f_array((amp_f_array>=gamma(1)&amp_f_array<=gamma(2)))};
%alpha
amp_freqs(1,2) = {'hga'};
amp_freqs(2,2) = {hga_idx};
amp_freqs(3,2) = {amp_f_array((amp_f_array>hga(1)&amp_f_array<=hga(2)))};
%beta
amp_freqs(1,3) = {'broadgamma'};
amp_freqs(2,3) = {broadgamma_idx};
amp_freqs(3,3) = {amp_f_array((amp_f_array>=broadgamma(1)&amp_f_array<=broadgamma(2)))};


%%%%extract plv values by frequency band
    %for each phase freq of interest - need to extract plvs for size phase x amp
    %freq plvs (for each phase freq of interest should define 3 blocks)
    
% delta_theta phase chunks 
    % means by subj
delta_theta_gamma_means = [];
delta_theta_hga_means = [];
delta_theta_broadgamma_means = [];
    % block plv mat by subj
delta_theta_gamma_mat = {};
delta_theta_hga_mat = {};
delta_theta_broadgamma_mat = {};
    %sig only
sig_delta_theta_gamma_means = [];
sig_delta_theta_hga_means = [];
sig_delta_theta_broadgamma_means = [];
sig_delta_theta_gamma_mat = {};
sig_delta_theta_hga_mat = {};
sig_delta_theta_broadgamma_mat = {};
    %nonsig only
nonsig_delta_theta_gamma_means = [];
nonsig_delta_theta_hga_means = [];
nonsig_delta_theta_broadgamma_means = [];
nonsig_delta_theta_gamma_mat = {};
nonsig_delta_theta_hga_mat = {};
nonsig_delta_theta_broadgamma_mat = {};

% alpha phase chunks 
    % means by subj
alpha_gamma_means = [];
alpha_hga_means = [];
alpha_broadgamma_means = [];
    % block plv mat by subj
alpha_gamma_mat = {};
alpha_hga_mat = {};
alpha_broadgamma_mat = {};
    %sig only
sig_alpha_gamma_means = [];
sig_alpha_hga_means = [];
sig_alpha_broadgamma_means = [];
sig_alpha_gamma_mat = {};
sig_alpha_hga_mat = {};
sig_alpha_broadgamma_mat = {};
    %nonsig only
nonsig_alpha_gamma_means = [];
nonsig_alpha_hga_means = [];
nonsig_alpha_broadgamma_means = [];
nonsig_alpha_gamma_mat = {};
nonsig_alpha_hga_mat = {};
nonsig_alpha_broadgamma_mat = {};

% beta phase chunks 
    % means by subj
beta_gamma_means = [];
beta_hga_means = [];
beta_broadgamma_means = [];
    % block plv mat by subj
beta_gamma_mat = {};
beta_hga_mat = {};
beta_broadgamma_mat = {};
    %sig only
sig_beta_gamma_means = [];
sig_beta_hga_means = [];
sig_beta_broadgamma_means = [];
sig_beta_gamma_mat = {};
sig_beta_hga_mat = {};
sig_beta_broadgamma_mat = {};
    %nonsig only
nonsig_beta_gamma_means = [];
nonsig_beta_hga_means = [];
nonsig_beta_broadgamma_means = [];
nonsig_beta_gamma_mat = {};
nonsig_beta_hga_mat = {};
nonsig_beta_broadgamma_mat = {};

% mean plv by subj - total, sig, non sig
all_elec_mean_plv_norm = [];
sig_elec_mean_plv_norm = [];
nonsig_elec_mean_plv_norm = [];


for s=1:length(subjs) %iterate through all subjects   
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    subj_sigs = data.sig_elecs; %getting labels for subject's sig elecs
    
    %extract mean all elec info 
    all_elec_mean_plv_norm(s) = data.mean_norm_plv_all_elec;
    
    %means for total data
    subj_delta_theta_gamma_mean = [];
    subj_delta_theta_hga_mean = [];
    subj_delta_theta_broadgamma_mean = [];
    
    subj_alpha_gamma_mean = [];
    subj_alpha_hga_mean = [];
    subj_alpha_broadgamma_mean = [];
    
    subj_beta_gamma_mean = [];
    subj_beta_hga_mean = [];
    subj_beta_broadgamma_mean = [];
    
    %means for sig only data
    sig_subj_delta_theta_gamma_mean = [];
    sig_subj_delta_theta_hga_mean = [];
    sig_subj_delta_theta_broadgamma_mean = [];
    
    sig_subj_alpha_gamma_mean = [];
    sig_subj_alpha_hga_mean = [];
    sig_subj_alpha_broadgamma_mean = [];
    
    sig_subj_beta_gamma_mean = [];
    sig_subj_beta_hga_mean = [];
    sig_subj_beta_broadgamma_mean = [];
    
    %means for nonsig data
    
    nonsig_subj_delta_theta_gamma_mean = [];
    nonsig_subj_delta_theta_hga_mean = [];
    nonsig_subj_delta_theta_broadgamma_mean = [];
    
    nonsig_subj_alpha_gamma_mean = [];
    nonsig_subj_alpha_hga_mean = [];
    nonsig_subj_alpha_broadgamma_mean = [];
    
    nonsig_subj_beta_gamma_mean = [];
    nonsig_subj_beta_hga_mean = [];
    nonsig_subj_beta_broadgamma_mean = [];

    se = 0; %count of sig elecs within subj
    nse = 0; %count of nonsig elec within subj
    for e=1:data.n_elecs
        
        elec = data.ofc_elecs(e);
        plv = data.norm_plvs{e}; %electrode plv is normalized plv matrix 

        %%%%extract plv values by phase frequency band
           %extract for total elecs
        %delta_theta
        delta_theta_gamma_plv = plv(gamma_idx,delta_theta_idx); %cluster matrix
        delta_theta_gamma_mat{1,s}(1,e) = {delta_theta_gamma_plv}; %add elec matrix to all subj array
        delta_theta_gamma_mean = mean(delta_theta_gamma_plv,'all'); %mean plv of cluster 
        subj_delta_theta_gamma_mean(e) = delta_theta_gamma_mean; %add elec mean to subj mean vector
        
        delta_theta_hga_plv = plv(hga_idx,delta_theta_idx); %cluster matrix
        delta_theta_hga_mat{1,s}(1,e) = {delta_theta_hga_plv}; %add elec matrix to all subj array
        delta_theta_hga_mean = mean(delta_theta_hga_plv,'all'); %mean plv of cluster 
        subj_delta_theta_hga_mean(e) = delta_theta_hga_mean; %add elec mean to subj mean vector
        
        delta_theta_broadgamma_plv = plv(broadgamma_idx,delta_theta_idx);
        delta_theta_broadgamma_mat{1,s}(1,e) = {delta_theta_broadgamma_plv};
        delta_theta_broadgamma_mean = mean(delta_theta_broadgamma_plv,'all');
        subj_delta_theta_broadgamma_mean(e) = delta_theta_broadgamma_mean;

        %alpha
        alpha_gamma_plv = plv(gamma_idx,alpha_idx);
        alpha_gamma_mat{1,s}(1,e) = {alpha_gamma_plv};
        alpha_gamma_mean = mean(alpha_gamma_plv,'all'); %mean plv of total cluster 
        subj_alpha_gamma_mean(e) = alpha_gamma_mean;
        
        alpha_hga_plv = plv(hga_idx,alpha_idx);
        alpha_hga_mat{1,s}(1,e) = {alpha_hga_plv};
        alpha_hga_mean = mean(alpha_hga_plv,'all'); %mean plv of total cluster 
        subj_alpha_hga_mean(e) = alpha_hga_mean;
        
        alpha_broadgamma_plv = plv(broadgamma_idx,alpha_idx);
        alpha_broadgamma_mat{1,s}(1,e) = {alpha_broadgamma_plv};
        alpha_broadgamma_mean = mean(alpha_broadgamma_plv,'all'); %mean plv of total cluster 
        subj_alpha_broadgamma_mean(e) = alpha_broadgamma_mean;

        %beta
        beta_gamma_plv = plv(gamma_idx,beta_idx);
        beta_gamma_mat{1,s}(1,e) = {beta_gamma_plv};
        beta_gamma_mean = mean(beta_gamma_plv,'all'); %mean plv of total cluster 
        subj_beta_gamma_mean(e) = beta_gamma_mean; %add to subj mean vector
        
        beta_hga_plv = plv(hga_idx,beta_idx);
        beta_hga_mat{1,s}(1,e) = {beta_hga_plv};
        beta_hga_mean = mean(beta_hga_plv,'all');
        subj_beta_hga_mean(e) = beta_hga_mean;
        
        beta_broadgamma_plv = plv(broadgamma_idx,beta_idx);
        beta_broadgamma_mat{1,s}(1,e) = {beta_broadgamma_plv};
        beta_broadgamma_mean = mean(beta_broadgamma_plv,'all'); %mean plv of total cluster 
        subj_beta_broadgamma_mean(e) = beta_broadgamma_mean;
        
        % add info to either sig or nonsig data
        % check whether electrode is significant
        is_sig = sum(find(subj_sigs==elec)~=0); %check if elec is in sig list - sum of nonzero values should be nonzero if elec is sig, if nonsig sum =0
        % SIG = is_sig ~= 0 ; NONSIG = is_sig == 0
        
        if is_sig ~= 0 %significant electrode
            se=se+1; %count of sig elecs within subj
            
            sig_elec_mean_plv_norm(se) = data.sig_mean_norm_plv;
            
            %delta-theta
            sig_subj_delta_theta_gamma_mean(se) = delta_theta_gamma_mean;
            sig_subj_delta_theta_hga_mean(se) = delta_theta_hga_mean;
            sig_subj_delta_theta_broadgamma_mean(se) = delta_theta_broadgamma_mean;
            sig_delta_theta_gamma_mat{1,s}(1,se) = {delta_theta_gamma_plv};
            sig_delta_theta_hga_mat{1,s}(1,se) = {delta_theta_hga_plv};
            sig_delta_theta_broadgamma_mat{1,s}(1,se) = {delta_theta_broadgamma_plv};
            %alpha    
            sig_subj_alpha_gamma_mean(se) = alpha_gamma_mean;
            sig_subj_alpha_hga_mean(se) = alpha_hga_mean;
            sig_subj_alpha_broadgamma_mean(se) = alpha_broadgamma_mean;
            sig_alpha_gamma_mat{1,s}(1,se) = {alpha_gamma_plv};
            sig_alpha_hga_mat{1,s}(1,se) = {alpha_hga_plv};
            sig_alpha_broadgamma_mat{1,s}(1,se) = {alpha_broadgamma_plv};
            %beta
            sig_subj_beta_gamma_mean(se) = beta_gamma_mean;
            sig_subj_beta_hga_mean(se) = beta_hga_mean;
            sig_subj_beta_broadgamma_mean(se) = beta_broadgamma_mean;
            sig_beta_gamma_mat{1,s}(1,se) = {beta_gamma_plv};
            sig_beta_hga_mat{1,s}(1,se) = {beta_hga_plv};
            sig_beta_broadgamma_mat{1,s}(1,se) = {beta_broadgamma_plv};
        
        elseif is_sig == 0 %nonsignificant electrode 
            nse=nse+1; %count of nonsig elecs within subj
            
            nonsig_elec_mean_plv_norm(nse) = data.nonsig_mean_norm_plv;
            
            %delta-theta
            nonsig_subj_delta_theta_gamma_mean(nse) = delta_theta_gamma_mean;
            nonsig_subj_delta_theta_hga_mean(nse) = delta_theta_hga_mean;
            nonsig_subj_delta_theta_broadgamma_mean(nse) = delta_theta_broadgamma_mean;
            nonsig_delta_theta_gamma_mat{1,s}(1,nse) = {delta_theta_gamma_plv};
            nonsig_delta_theta_hga_mat{1,s}(1,nse) = {delta_theta_hga_plv};
            nonsig_delta_theta_broadgamma_mat{1,s}(1,nse) = {delta_theta_broadgamma_plv};
            %alpha
            nonsig_subj_alpha_gamma_mean(nse) = alpha_gamma_mean;
            nonsig_subj_alpha_hga_mean(nse) = alpha_hga_mean;
            nonsig_subj_alpha_broadgamma_mean(nse) = alpha_broadgamma_mean;
            nonsig_alpha_gamma_mat{1,s}(1,nse) = {alpha_gamma_plv};
            nonsig_alpha_hga_mat{1,s}(1,nse) = {alpha_hga_plv};
            nonsig_alpha_broadgamma_mat{1,s}(1,nse) = {alpha_broadgamma_plv};
            %beta
            nonsig_subj_beta_gamma_mean(nse) = beta_gamma_mean;
            nonsig_subj_beta_hga_mean(nse) = beta_hga_mean;
            nonsig_subj_beta_broadgamma_mean(nse) = beta_broadgamma_mean;
            nonsig_beta_gamma_mat{1,s}(1,nse) = {beta_gamma_plv};
            nonsig_beta_hga_mat{1,s}(1,nse) = {beta_hga_plv};
            nonsig_beta_broadgamma_mat{1,s}(1,nse) = {beta_broadgamma_plv};
            
        end 
        
        

    end %elec loop
    
    %mean across all elec per subj 
    delta_theta_gamma_means(1,s) = mean(subj_delta_theta_gamma_mean);
    delta_theta_hga_means(1,s) = mean(subj_delta_theta_hga_mean);
    delta_theta_broadgamma_means(1,s) = mean(subj_delta_theta_broadgamma_mean);
    
    alpha_gamma_means(1,s) = mean(subj_alpha_gamma_mean);
    alpha_hga_means(1,s) = mean(subj_alpha_hga_mean);
    alpha_broadgamma_means(1,s) = mean(subj_alpha_broadgamma_mean);
    
    beta_gamma_means(1,s) = mean(subj_beta_gamma_mean);
    beta_hga_means(1,s) = mean(subj_beta_hga_mean);
    beta_broadgamma_means(1,s) = mean(subj_beta_broadgamma_mean);
    
    %mean across SIG elec per subj 
    sig_delta_theta_gamma_means(1,s) = mean(sig_subj_delta_theta_gamma_mean);
    sig_delta_theta_hga_means(1,s) = mean(sig_subj_delta_theta_hga_mean);
    sig_delta_theta_broadgamma_means(1,s) = mean(sig_subj_delta_theta_broadgamma_mean);
    
    sig_alpha_gamma_means(1,s) = mean(sig_subj_alpha_gamma_mean);
    sig_alpha_hga_means(1,s) = mean(sig_subj_alpha_hga_mean);
    sig_alpha_broadgamma_means(1,s) = mean(sig_subj_alpha_broadgamma_mean);
    
    sig_beta_gamma_means(1,s) = mean(sig_subj_beta_gamma_mean);
    sig_beta_hga_means(1,s) = mean(sig_subj_beta_hga_mean);
    sig_beta_broadgamma_means(1,s) = mean(sig_subj_beta_broadgamma_mean);
    
    %mean across NONSIG elec per subj 
    nonsig_delta_theta_gamma_means(1,s) = mean(nonsig_subj_delta_theta_gamma_mean);
    nonsig_delta_theta_hga_means(1,s) = mean(nonsig_subj_delta_theta_hga_mean);
    nonsig_delta_theta_broadgamma_means(1,s) = mean(nonsig_subj_delta_theta_broadgamma_mean);
    
    nonsig_alpha_gamma_means(1,s) = mean(nonsig_subj_alpha_gamma_mean);
    nonsig_alpha_hga_means(1,s) = mean(nonsig_subj_alpha_hga_mean);
    nonsig_alpha_broadgamma_means(1,s) = mean(nonsig_subj_alpha_broadgamma_mean);
    
    nonsig_beta_gamma_means(1,s) = mean(nonsig_subj_beta_gamma_mean);
    nonsig_beta_hga_means(1,s) = mean(nonsig_subj_beta_hga_mean);
    nonsig_beta_broadgamma_means(1,s) = mean(nonsig_subj_beta_broadgamma_mean);
        
        
end %subj loop
    %%%% populate all elec structure
    %add cluster indexing info to all_elec_plv_info
    all_elec_plv_info.subjs = subjs;
    all_elec_plv_info.mean_plv_norm = all_elec_mean_plv_norm;
    
    %add delta-theta info
    all_elec_plv_info.delta_theta_gamma_means = delta_theta_gamma_means;
    all_elec_plv_info.delta_theta_hga_means = delta_theta_hga_means;
    all_elec_plv_info.delta_theta_broadgamma_means = delta_theta_broadgamma_means;
    all_elec_plv_info.delta_theta_gamma_plvs = delta_theta_gamma_mat;
    all_elec_plv_info.delta_theta_hga_plvs = delta_theta_hga_mat;
    all_elec_plv_info.delta_theta_broadgamma_plvs = delta_theta_broadgamma_mat;
    
    %add alpha info
    all_elec_plv_info.alpha_gamma_means = alpha_gamma_means;
    all_elec_plv_info.alpha_hga_means = alpha_hga_means;
    all_elec_plv_info.alpha_broadgamma_means = alpha_broadgamma_means;
    all_elec_plv_info.alpha_gamma_plvs = alpha_gamma_mat;
    all_elec_plv_info.alpha_hga_plvs = alpha_hga_mat;
    all_elec_plv_info.alpha_broadgamma_plvs = alpha_broadgamma_mat;
    
    %add beta info
    all_elec_plv_info.beta_gamma_means = beta_gamma_means;
    all_elec_plv_info.beta_hga_means = beta_hga_means;
    all_elec_plv_info.beta_broadgamma_means = beta_broadgamma_means;
    all_elec_plv_info.beta_gamma_plvs = beta_gamma_mat;
    all_elec_plv_info.beta_hga_plvs = beta_hga_mat;
    all_elec_plv_info.beta_broadgamma_plvs = beta_broadgamma_mat;
    
    %%%% populate sig elec structure
    %add cluster indexing info to SIG_elec_plv_info
    sig_elec_plv_info.mean_plv_norm = sig_elec_mean_plv_norm;
    
    %add delta-theta info
    sig_elec_plv_info.delta_theta_gamma_means = sig_delta_theta_gamma_means;
    sig_elec_plv_info.delta_theta_hga_means = sig_delta_theta_hga_means;
    sig_elec_plv_info.delta_theta_broadgamma_means = sig_delta_theta_broadgamma_means;
    sig_elec_plv_info.delta_theta_gamma_plvs = sig_delta_theta_gamma_mat;
    sig_elec_plv_info.delta_theta_hga_plvs = sig_delta_theta_hga_mat;
    sig_elec_plv_info.delta_theta_broadgamma_plvs = sig_delta_theta_broadgamma_mat;
    
    %add alpha info
    sig_elec_plv_info.alpha_gamma_means = sig_alpha_gamma_means;
    sig_elec_plv_info.alpha_hga_means = sig_alpha_hga_means;
    sig_elec_plv_info.alpha_broadgamma_means = sig_alpha_broadgamma_means;
    sig_elec_plv_info.alpha_gamma_plvs = sig_alpha_gamma_mat;
    sig_elec_plv_info.alpha_hga_plvs = sig_alpha_hga_mat;
    sig_elec_plv_info.alpha_broadgamma_plvs = sig_alpha_broadgamma_mat;
    
    %add beta info
    sig_elec_plv_info.beta_gamma_means = sig_beta_gamma_means;
    sig_elec_plv_info.beta_hga_means = sig_beta_hga_means;
    sig_elec_plv_info.beta_broadgamma_means = sig_beta_broadgamma_means;
    sig_elec_plv_info.beta_gamma_plvs = sig_beta_gamma_mat;
    sig_elec_plv_info.beta_hga_plvs = sig_beta_hga_mat;
    sig_elec_plv_info.beta_broadgamma_plvs = sig_beta_broadgamma_mat;
    
    %%%% populate NONsig elec structure
    %add cluster indexing info to all_elec_plv_info
    nonsig_elec_plv_info.mean_plv_norm = nonsig_elec_mean_plv_norm;
    
    %add delta-theta info
    nonsig_elec_plv_info.delta_theta_gamma_means = nonsig_delta_theta_gamma_means;
    nonsig_elec_plv_info.delta_theta_hga_means = nonsig_delta_theta_hga_means;
    nonsig_elec_plv_info.delta_theta_broadgamma_means = nonsig_delta_theta_broadgamma_means;
    nonsig_elec_plv_info.delta_theta_gamma_plvs = nonsig_delta_theta_gamma_mat;
    nonsig_elec_plv_info.delta_theta_hga_plvs = nonsig_delta_theta_hga_mat;
    nonsig_elec_plv_info.delta_theta_broadgamma_plvs = nonsig_delta_theta_broadgamma_mat;
    
    %add alpha info
    nonsig_elec_plv_info.alpha_gamma_means = nonsig_alpha_gamma_means;
    nonsig_elec_plv_info.alpha_hga_means = nonsig_alpha_hga_means;
    nonsig_elec_plv_info.alpha_broadgamma_means = nonsig_alpha_broadgamma_means;
    nonsig_elec_plv_info.alpha_gamma_plvs = nonsig_alpha_gamma_mat;
    nonsig_elec_plv_info.alpha_hga_plvs = nonsig_alpha_hga_mat;
    nonsig_elec_plv_info.alpha_broadgamma_plvs = nonsig_alpha_broadgamma_mat;
    
    %add beta info
    nonsig_elec_plv_info.beta_gamma_means = nonsig_beta_gamma_means;
    nonsig_elec_plv_info.beta_hga_means = nonsig_beta_hga_means;
    nonsig_elec_plv_info.beta_broadgamma_means = nonsig_beta_broadgamma_means;
    nonsig_elec_plv_info.beta_gamma_plvs = nonsig_beta_gamma_mat;
    nonsig_elec_plv_info.beta_hga_plvs = nonsig_beta_hga_mat;
    nonsig_elec_plv_info.beta_broadgamma_plvs = nonsig_beta_broadgamma_mat;


end %functionend
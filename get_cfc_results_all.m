function [cfc_all_subj] = get_cfc_results_all(subjs,data_path)

%inputs: subj id, data path to plv matrix & surrogate plv matrix
%output: all subj structure containing raw plv matrix, surrogate plv
    %matrix, normalized plv matrix, mean plvs by and across electrodes,
    %pixel-level pvalue matrices (proportion and zstatistic)

for s=1:length(subjs)   
    
    subj_id = char(subjs(s));

    
    elec = load(strcat(data_path,subj_id,'/',subj_id,'_elec.mat'));
    anat_info = elec.elec.anat;
    ofc_idx = find(ismember(anat_info,'OrG'));
    n_ofc_elecs = length(ofc_idx);

    subj_results_path = [data_path '/cfc_results/' subj_id '/'];
    plv_str = [subj_id '_plv_matrix_real_'];
    perm_str = [subj_id '_plv_matrix_perm_'];
    
    %cell arrays for full matrices/elec
    real_plvs = {};
    perm_plvs = {};
    norm_plvs = {};
    %cell arrays for mean values/elec
    mean_real_plv_by_elec = {};
    mean_norm_plv_by_elec = {};
    %matrices for mean values/subj (all elecs)
    mean_real_plv_all_elec = []; %must mean all values before structure assignment
    mean_norm_plv_all_elec = []; %must mean all values before structure assignment
    %cell arrays for statistics
    prop_pvals_by_pixel = {}; %pvalue = (prop surr > pixel_plv)/num_surr
    zstat_pvals_by_pixel = {}; %pvalue from zscore of pixel = (1-normcdf(zscore_pixel));
    
    for e=1:length(ofc_idx)
        
        elec_num = ofc_idx(e); 
        load(strcat(subj_results_path,plv_str,tostring(elec_num),'.mat'))
        load(strcat(subj_results_path,perm_str,tostring(elec_num),'.mat'))

        
        real_plvs{1,e} = plv_matrix_real;
        perm_plvs{1,e} = plv_matrix_perm;   
        
        mean_real_plv_by_elec{1,e} = mean(plv_matrix_real,'all'); %mean of raw plv for electrode 
        mean_real_plv_all_elec(1,e) = mean(plv_matrix_real,'all'); %mean of raw plv for electrode > will be group meaned
        
        
        elec_norm_plv = zeros(size(plv_matrix_real)); %this needs to be within pixel based on permutation data!
        elec_prop_pvals = zeros(size(plv_matrix_real));
        elec_zstat_pvals = zeros(size(plv_matrix_real));

        
        %iterating through pixels in each df 
        for a=1:40 %row-wise amplitudes
            for p=1:19 %column-wise phases
                
                %extract pixel data
                real_pixel = plv_matrix_real(a,p); %real plv for pixel
                perm_pixels = squeeze(plv_matrix_perm(a,p,:)); %surrogate data for pixel
                num_perms =  length(perm_pixels); %number of permutations
                
                %normalize pixel data
                perm_mean = mean(perm_pixels); %mean of pixel permutation distribution
                perm_std = std(perm_pixels); %std of pixel permutation distribution
                zscore_pixel = (real_pixel-perm_mean)/perm_std; %normalized pixel value from permutation distribution
                elec_norm_plv(a,p) = zscore_pixel; %add normalized pixel value to within elec norm plv matrix
                
                %%%statistical tests%%%
                
                %pvalue from proportion (pn):
                    %calculate proportion of surrogate pixels > real pixel value NOT ZSCORED PIXEL!
                n_perm_greater = sum(perm_pixels>real_pixel); %num surrogates greater than real pixel plv
                elec_prop_pvals(a,p) = (n_perm_greater/num_perms); %pval calculation as proportion
                
                %pvalue from zstatistic(pz):
                    %pixel zscore is zstatistic for normcdf
                    %permutation dist, calculate pvalue from zscore
                zstat_pval = (1-normcdf(zscore_pixel));
                elec_zstat_pvals(a,p) = zstat_pval; %pvalue calculation with t statistic, n-1 df
            
            
            end 
            

        end %end of single elec all pixel iteration
        
        %store electrode norm plv mat information 
        norm_plvs{1,e} = elec_norm_plv;
        mean_norm_plv_by_elec{1,e} = mean(elec_norm_plv,'all'); %store mean norm plv for each elec
        mean_norm_plv_all_elec(1,e) = mean(elec_norm_plv,'all'); %to be averaged across elec by subj
        
        %store electrode stats (uncorrected) 
        prop_pvals_by_pixel{1,e} = elec_prop_pvals; %pvalue = (prop surr > pixel_plv)/num_surr
        zstat_pvals_by_pixel{1,e} = elec_zstat_pvals;%pvalue from tstatistic by pixel (tcdf(t_stat,(length(perm_pixels)-1),'upper'));


    end %end of all elec iteration
    
    %store single subj all elec information in subj_struct
    subj_struct.subj_id = subj_id; %subj_id
    subj_struct.ofc_elecs = ofc_idx; %ofc electrode ids
    subj_struct.n_elecs = n_ofc_elecs; %number of ofc elecs
    subj_struct.real_plvs = real_plvs; %all raw plvs by elec
    subj_struct.perm_plvs = perm_plvs; %all permutation distributions by elec
    subj_struct.norm_plvs = norm_plvs; %all normalized plvs by elec
    subj_struct.mean_real_plv_by_elec = mean_real_plv_by_elec; %mean raw plv by electrode
    subj_struct.mean_real_plv_all_elec = mean(mean_real_plv_all_elec,'all'); %mean raw plv for ALL electrodes (1 value/subj)
    subj_struct.mean_norm_plv_by_elec = mean_norm_plv_by_elec; %mean of norm plv by electrode
    subj_struct.mean_norm_plv_all_elec = mean(mean_norm_plv_all_elec,'all'); %mean of norm plv for ALL electrodes (1 value/subj)
    subj_struct.prop_pvals = prop_pvals_by_pixel; %matrix of proportion pvalues for each pixel by electrode
    subj_struct.zstat_pvals = zstat_pvals_by_pixel;%matrix of ttest pvalues for each pixel by electrode

    %assign subj_struct to all subj struct
    cfc_all_subj.(subj_id) = subj_struct;
    

end %end of subj iteration



end 

%last updated 10/26/2023 - changed p value calculations to reflect ANTS
%cohen pn and pz 


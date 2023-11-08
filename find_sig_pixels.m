function [sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count,...
    find_sig_elecs_header,cfc_all_subj] = find_sig_pixels(cfc_all_subj,subjs,pval_type)



% header with info for sig elecs and non sig elecs data 
find_sig_elecs_header.pval_type = {pval_type}; %which significance test generated pvals (pval_type either 'prop' or 'zstat') 
find_sig_elecs_header.row_labels = {'num_elecs','elec_ids','all_elec_mean_raw_plv','all_elec_mean_norm_plv',...
    'adj_pvals','adj_norm_plv','num_sig_pixels'}; %nonsig_elecs_by_subj don't have last two 


sig_elecs_by_subj = {}; %row corresponds to find_sig_elecs_header.row_labels
nonsig_elecs_by_subj = {}; %row data corresponds to find_sig_elecs_header.row_labels


sig_count = 0; %counter for total number of sig electrodes across subj
nonsig_count = 0; %counter for total number of nonsig electrodes across subj

for s=1:length(subjs)   
    
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    
    subj_count=0; %num significant electrodes within subj
    subj_sig_elecs = []; %elec ids of subj sig electrodes
    subj_sig_mean_raw_plv = []; % raw plv means for all sig elecs - will be averaged
    subj_sig_mean_norm_plv = []; % norm plv means for all sig elecs - will be averaged
    subj_sig_adj_pvals = {}; %adjusted pvals from fdr correction for each electrode
    subj_sig_adj_norm_plv = {}; %adjusted plv mat only sig pixels from fdr correction for each electrode
    subj_num_sig_pixel = {}; %num sig pixels for each electrode
    
    nonsig_subj_count=0; %num significant electrodes within subj
    subj_nonsig_elecs = []; %elec ids of subj nonsig electrodes
    subj_nonsig_mean_raw_plv = []; % raw plv means for all nonsig elecs - will be averaged
    subj_nonsig_mean_norm_plv = []; % norm plv means for all nonsig elecs - will be averaged
    subj_nonsig_adj_pvals = {}; %adjusted pvals from fdr correction for each nonsig electrode
    
    for e=1:data.n_elecs
        
                
        %%%% To filter PLV mats to include only significant coupling values
        %step 1: extract pixel pvalues for each PLV matrix 
        %step 2: FDR correction on p values 
        %step 3: set all pixels with p value >= 0.05 to zero
        %step 4: if there are no significant pixels > elec is non
            %significant
        %step 5: save corrected normalized PLV matrix into cfc all subj
        
        %extract electrode data 
        elec = data.ofc_elecs(e);
        elec_mean_raw_plv = data.mean_real_plv_by_elec{e}; %mean raw plv for elec
        elec_mean_norm_plv = data.mean_norm_plv_by_elec{e}; %mean norm plv for elec 
        plv = data.norm_plvs{e}; %electrode plv is normalized plv matrix 
        
        %step 1:
        %get pvalues for fdr correction - depends on pval_type input
        if strcmp(pval_type,'prop') 
            raw_pvals = data.prop_pvals{e}; %pvalues calculated from perm for each pixel (#perm > pixel/# perms total)
        elseif strcmp(pval_type, 'zstat')
            raw_pvals = data.zstat_pvals{e}; %pvalues calculated from perm for each pixel (t test)
        end 
        
        %step 2: 
        
        %Benjamani & Hochberg False Discovery Rate (FDR) Correction (function = fdr_bh, 2015, David Groppe)
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(raw_pvals,0.05,'dep'); %extract corrected pvalue matrix from fdr correction

        
        %step 3:
        [row,col] = find(adj_p < 0.05); %find indices of sig plvs (corrected pval < 0.05)
        sig_plv = plv(row,col); %extract norm plv values for sig pixels
        adj_plv = zeros(size(plv)); %make empty adj plv matrix
        adj_plv(row,col) = sig_plv; %include plvs for only sig pixels
        
        %step 4:
        if sum(adj_plv,'all') == 0 %if there are no significant pixels 
            % add data to nonsig_elecs_by_subj
            nonsig_subj_count = nonsig_subj_count+1; %add one to num nonsig elecs within subj
            nonsig_count = nonsig_count+1; %add one to total num nonsig elecs
            subj_nonsig_elecs(nonsig_subj_count,:) = elec; %add elec id to nonsig elec tracker
            subj_nonsig_mean_raw_plv(nonsig_subj_count,:) = elec_mean_raw_plv; % raw plv mean
            subj_nonsig_mean_norm_plv(nonsig_subj_count,:) = elec_mean_norm_plv; % norm plv mean
            subj_nonsig_adj_pvals(nonsig_subj_count,:) = {adj_p}; %adjusted pvals from fdr correction
            
 
        else % add data to sig_elecs_by_subj
            subj_count=subj_count+1; %add one to num sig elecs within subj
            sig_count = sig_count+1; %add one to total num sig elecs
            subj_sig_elecs(subj_count,:) = elec; %add elec id to sig elec tracker
            subj_sig_mean_raw_plv(subj_count,:) = elec_mean_raw_plv; % raw plv mean
            subj_sig_mean_norm_plv(subj_count,:) = elec_mean_norm_plv; % norm plv mean
            subj_sig_adj_pvals(subj_count,:) = {adj_p}; %adjusted pvals from fdr correction
            subj_sig_adj_norm_plv(subj_count,:) = {adj_plv}; % corrected plv mat - sig pixels only!
            subj_num_sig_pixel(subj_count,:) = {sum(adj_plv~=0,'all')}; % num sig pixels!

        end 
        
    end %end of electrode loop
    
    %calculate means of all sig/nonsig plv - must calculate b4 assignment
    sig_mean_raw_plv =  mean(subj_sig_mean_raw_plv,'all'); %all_elec_mean_raw_plv
    sig_mean_norm_plv = mean(subj_sig_mean_norm_plv,'all'); %all_elec_mean_norm_plv
    nonsig_mean_raw_plv = mean(subj_nonsig_mean_raw_plv,'all'); %all_elec_mean_raw_plv
    nonsig_mean_norm_plv = mean(subj_nonsig_mean_norm_plv,'all'); %all_elec_mean_norm_plv
    
    %within subj summary
    sig_elecs_by_subj(1,s) = {subj_count}; %num sig elecs within subj
    sig_elecs_by_subj(2,s) = {subj_sig_elecs}; %elec id of sig elecs within subj
    sig_elecs_by_subj(3,s) = {sig_mean_raw_plv}; %all_elec_mean_raw_plv
    sig_elecs_by_subj(4,s) = {sig_mean_norm_plv}; %all_elec_mean_norm_plv
    sig_elecs_by_subj(5,s) = {subj_sig_adj_pvals}; %adjusted pvals from fdr correction for each sig elec
    sig_elecs_by_subj(6,s) = {subj_sig_adj_norm_plv}; % corrected plv mats
    sig_elecs_by_subj(7,s) = {subj_num_sig_pixel}; % num sig pixels


    nonsig_elecs_by_subj(1,s) = {nonsig_subj_count}; %num nonsig elecs within subj
    nonsig_elecs_by_subj(2,s) = {subj_nonsig_elecs}; %elec id of nonsig elecs within subj
    nonsig_elecs_by_subj(3,s) = {nonsig_mean_raw_plv}; %all_elec_mean_raw_plv
    nonsig_elecs_by_subj(4,s) = {nonsig_mean_norm_plv}; %all_elec_mean_norm_plv
    nonsig_elecs_by_subj(5,s) = {subj_nonsig_adj_pvals}; %adjusted pvals from fdr correction for each nonsig elec

    %data to update all subj struct  
        %add sig data
    data.sig_elec_num = subj_count; %num sig elecs for each subj
    data.sig_elecs = subj_sig_elecs; % elec ids of sig elecs
    data.sig_mean_real_plv = sig_mean_raw_plv; %mean raw plv of sig elecs
    data.sig_mean_norm_plv = sig_mean_norm_plv; %mean norm plv of sig elecs
    
        %add non significant data
    data.nonsig_elec_num = nonsig_subj_count; %num nonsig elecs for each subj
    data.nonsig_elecs = subj_nonsig_elecs; % elec ids of nonsig elecs
    data.nonsig_mean_real_plv = nonsig_mean_raw_plv; %mean raw plv of nonsig elecs
    data.nonsig_mean_norm_plv = nonsig_mean_norm_plv; %mean norm plv of nonsig elecs

    % update cfc_all_subj subj field updated data structure
    cfc_all_subj.(subj_id) = data;
    

    end %end of all subj loop




function [sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count,...
    find_sig_elecs_header,cfc_all_subj] = find_sig_elecs(cfc_all_subj,subjs,pval_type,plot_elec)

%pval_type can either be 'prop' or 'zstat'
%plot_elec can either be 'yes' or 'no'

% header with info for sig elecs and non sig elecs data 
find_sig_elecs_header.pval_type = {pval_type}; %which significance test generated pvals (pval_type either 'prop' or 'zstat') 
find_sig_elecs_header.row_labels = {'num_elecs','elec_ids','all_elec_mean_raw_plv','all_elec_mean_norm_plv',...
    'elec_norm_plv_mats','adj_pvals','cluster_amps','cluster_phases','cluster_mat','cluster_pvals','cluster_sig_pixel_prop'};


sig_elecs_by_subj = {}; %row corresponds to find_sig_elecs_header.row_labels
nonsig_elecs_by_subj = {}; %row data corresponds to find_sig_elecs_header.row_labels


sig_count = 0; %counter for total number of sig electrodes across subj
nonsig_count = 0; %counter for total number of nonsig electrodes across subj

for s=1:length(subjs)   
    
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    
    all_adj_pval_mats = {};
    
    subj_count=0; %num significant electrodes within subj
    subj_sig_elecs = []; %elec ids of subj sig electrodes
    subj_sig_mean_raw_plv = []; % raw plv means for all sig elecs - will be averaged
    subj_sig_mean_norm_plv = []; % norm plv means for all sig elecs - will be averaged
    subj_sig_norm_plv_mats = {}; % norm plv full matrices for all sig elecs
    subj_sig_adj_pvals = {}; %adjusted pvals from fdr correction for each electrode
    subj_sig_clust_amps = {}; %amplitude freqs of largest cluster for each electrode
    subj_sig_clust_phases = {}; %amplitude freqs of largest cluster for each electrode
    subj_sig_clust_mats = {}; %matrix of normalized plvs for most significant clusters
    subj_sig_clust_pvals = {}; %matrix of pvals corresponding to pixels in most significant clusters
    subj_sig_cluster_sig_pixel_prop = {}; %proportion of cluster pvalues <= 0.05 (for significant elecs should always >= 0.95)
    
    nonsig_subj_count=0; %num significant electrodes within subj
    subj_nonsig_elecs = []; %elec ids of subj nonsig electrodes
    subj_nonsig_mean_raw_plv = []; % raw plv means for all nonsig elecs - will be averaged
    subj_nonsig_mean_norm_plv = []; % norm plv means for all nonsig elecs - will be averaged
    subj_nonsig_norm_plv_mats = {}; % norm plv full matrices for all nonsig elecs
    subj_nonsig_adj_pvals = {}; %adjusted pvals from fdr correction for each nonsig electrode
    subj_nonsig_clust_amps = {}; %amplitude freqs of largest cluster for each nonsig electrode
    subj_nonsig_clust_phases = {}; %amplitude freqs of largest cluster for each nonsig electrode
    subj_nonsig_clust_mats = {}; %matrix of normalized plvs for largest cluster (nonsig)
    subj_nonsig_clust_pvals = {}; %matrix of pvals corresponding to pixels for largest cluster (nonsig)
    subj_nonsig_cluster_sig_pixel_prop = {}; %proportion of cluster pvalues <= 0.05 (if no cluster = 0)
    
    for e=1:data.n_elecs
        
        %extract electrode data 
        elec = data.ofc_elecs(e);
        elec_mean_raw_plv = data.mean_real_plv_by_elec{e}; %mean raw plv for elec
        elec_mean_norm_plv = data.mean_norm_plv_by_elec{e}; %mean norm plv for elec 
        plv = data.norm_plvs{e}; %electrode plv is normalized plv matrix 
        
        %get pvalues for fdr correction - depends on pval_type input
        if strcmp(pval_type,'prop') 
            raw_pvals = data.prop_pvals{e}; %pvalues calculated from perm for each pixel (#perm > pixel/# perms total)
        elseif strcmp(pval_type, 'zstat')
            raw_pvals = data.zstat_pvals{e}; %pvalues calculated from perm for each pixel (t test)
        end 
        %Benjamani & Hochberg False Discovery Rate (FDR) Correction (function = fdr_bh, 2015, David Groppe)
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(raw_pvals,0.05,'dep'); %extract corrected pvalue matrix from fdr correction

        %%%% To determine if an electrode has "significant cfc" 
            %step 1: extract contours from normalized plv matrix
            %step 2: determine cut off for significant clusters (depends on max zscore value)
            %step 3: extract phase/amp freq indices from sig clusters
            %step 4: correct freq indices to include 2x2 or larger clusters
            %step 5: extract adjusted pvalues for all pixels in cluster
            %step 6: determine if cluster is significant:
                        %proportion of corrected pixel pvals < 0.05 (alpha)
                        %must be greater than 0.95 (95% of clust pixels are significant)
        
        %step 1:       
        figure;
        [cm,h] = contourf(plv); %get contours of normalized plv mat
        if strcmp(plot_elec,'yes')
            single_elec_fig_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/data/cfc_data/cfc_results/';
            fig_save_path = [single_elec_fig_path subj_id '/'];
            amp_f_array = 5:5:200;
            phase_f_array = 2:20;
            % Set locations of ticks and their labels.
            set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
            'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',22,...
            'TickDir','in','FontName','Lucida Grande');
            title(['\rm' subj_id ' Electrode ' num2str(elec)],'FontSize',22,'FontName','Lucida Grande');
            xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Lucida Grande');
            ylabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Lucida Grande');
            h = colorbar;
            set(get(h, 'ylabel'), 'string', '\fontsize{26}PLV Z-Score','FontSize',22,'FontName','Lucida Grande');
            saveas(gca,[fig_save_path, subj_id, '_', num2str(elec), '_plv.pdf'])
        end 
        
        [contourTable,contourArray] = getContourLineCoordinates(cm); %extract contour information %function from 2020, Adam Danz
        %number of contours is max num Group in contourTable or len
        %contourArray
        
        %step 2:
        
        %option 1 - extract all clusters 2 std above mean 
        
%         max_z = 2; %2 std above mean (alpha 0.05)
%         max_level = contourTable.Level> max_z; %extract contours > 2 std above mean of perm dist
       
        %option 2 - extract clusters variably based on sig level
        if max(contourTable.Level) <= 2 %some elecs have a maximum zscore of 2 or lower
            max_z = 2; %2 std above mean (alpha 0.05)
            max_level = contourTable.Level >= max_z; %extract significant contour info
        else %(max(contourTable.Level) > 2 && max(contourTable.Level)<=5) %for clusters that have max std greater than 2 but less than/= 5
            max_z = max(contourTable.Level); %setting max to max cluster
            max_level = contourTable.Level >= max_z; %extracts all values of highest level clusters 
        end 

        %step 3:
        if sum(max_level) == 0 %some electrodes have no significant clusters (max level <2)
            %define null matrices for nonsig data structure
            cluster_mat = []; % no significant clusters
            cluster_pvals = []; % no pvalues - no cluster exists 
            cluster_sig_pixel_prop = 0; % no significant pixels
            sig_amp_idx = []; % no significant amps 
            sig_phase_idx = []; % no significant phases 
        
        else % electrodes that have contours > 2 std 
            
            %get indices of significant contours
            max_levelX = contourTable.X(max_level); %extract x index of max contours (phase freq)
            max_levelY = contourTable.Y(max_level); %extract y index of max contours (amp freq)
            sig_amp_idx = unique(round(max_levelY)); %indexes of contour amp freqs (rounded)
            sig_phase_idx = unique(round(max_levelX)); %indexes of contour phase freqs (rounded)
            %amp_f_array = 5:5:200; % values of amp freqs
            %phase_f_array = 2:20; % values of phase freqs
            %must round contour index to whole number (actual data indices can't be exact contour curve)

            %sig_amps = amp_f_array(sig_amp_idx'); %amp freqs corresponding to amp index
            %sig_phases = phase_f_array(sig_phase_idx'); %phase freqs corresponding to phase index
            
            %step 4:
                % to be considered a cluster - must have at least 2
                % consecutive significant frequency bands for phase and amp
                
            %uncomment to remove consecutive amplitude requirement   
%             cluster_mat = plv(sig_amp_idx,sig_phase_idx); %extract norm plv matrix of corrected cluster
%             cluster_pvals = adj_p(sig_amp_idx,sig_phase_idx); %get corrected pvalues of every pixel in cluster
%             cluster_sig_pixel_prop = sum(cluster_pvals<=0.05,'all')/numel(cluster_pvals); %calculate proportion of pvalues < 0.05 - if 95% of p values are < 0.05, cluster is significant!

            if sum(diff(sig_amp_idx)==1) == 0 %if there are no consecutive sig amps
                %define null matrices for nonsig data structure
                cluster_mat = []; % no significant clusters
                cluster_pvals = []; % no pvalues - no cluster exists 
                cluster_sig_pixel_prop = 0; % no significant pixels


            else 
%             %step 5: 
                cluster_mat = plv(sig_amp_idx,sig_phase_idx); %extract norm plv matrix of corrected cluster
                cluster_pvals = adj_p(sig_amp_idx,sig_phase_idx); %get corrected pvalues of every pixel in cluster
                cluster_sig_pixel_prop = sum(cluster_pvals<=0.05,'all')/numel(cluster_pvals); %calculate proportion of pvalues < 0.05 - if 95% of p values are < 0.05, cluster is significant!
                
            end 
                

        end
        

        if cluster_sig_pixel_prop >=0.95 % alpha = 0.05, if prop of sig cluster surpasses false positive rate - clust = significant 
            % add data to sig_elecs_by_subj
            
            subj_count=subj_count+1; %add one to num sig elecs within subj
            sig_count = sig_count+1; %add one to total num sig elecs
            subj_sig_elecs(subj_count,:) = elec; %add elec id to sig elec tracker
            subj_sig_mean_raw_plv(subj_count,:) = elec_mean_raw_plv; % raw plv mean
            subj_sig_mean_norm_plv(subj_count,:) = elec_mean_norm_plv; % norm plv mean
            subj_sig_norm_plv_mats(subj_count,:) = {plv}; % plv norm full matrix
            subj_sig_adj_pvals(subj_count,:) = {adj_p}; %adjusted pvals from fdr correction
            subj_sig_clust_amps(subj_count,:) = {sig_amp_idx}; % amp freq index (corrected!)
            subj_sig_clust_phases(subj_count,:) = {sig_phase_idx}; % phase freq index
            subj_sig_clust_mats(subj_count,:) = {cluster_mat}; %matrix of normalized plvs for most significant clusters
            subj_sig_clust_pvals(subj_count,:) = {cluster_pvals}; %matrix of pvals corresponding to pixels in most significant clusters
            subj_sig_cluster_sig_pixel_prop(subj_count,:) = {cluster_sig_pixel_prop}; %proportion of cluster pvalues <= 0.05
            all_adj_pval_mats(e,:) = {adj_p};
            
        else
            % add data to nonsig_elecs_by_subj
            nonsig_subj_count = nonsig_subj_count+1; %add one to num nonsig elecs within subj
            nonsig_count = nonsig_count+1; %add one to total num nonsig elecs
            subj_nonsig_elecs(nonsig_subj_count,:) = elec; %add elec id to nonsig elec tracker
            subj_nonsig_mean_raw_plv(nonsig_subj_count,:) = elec_mean_raw_plv; % raw plv mean
            subj_nonsig_mean_norm_plv(nonsig_subj_count,:) = elec_mean_norm_plv; % norm plv mean
            subj_nonsig_norm_plv_mats(nonsig_subj_count,:) = {plv}; % plv norm full matrix
            subj_nonsig_adj_pvals(nonsig_subj_count,:) = {adj_p}; %adjusted pvals from fdr correction
            subj_nonsig_clust_amps(nonsig_subj_count,:) = {sig_amp_idx}; %amp freqs of largest cluster (zero if no clusters)
            subj_nonsig_clust_phases(nonsig_subj_count,:) = {sig_phase_idx}; %phase freqs of largest cluster (zero if no clusters)
            subj_nonsig_clust_mats(nonsig_subj_count,:) = {cluster_mat}; %matrix of normalized plvs for most significant clusters (empty if none)
            subj_nonsig_clust_pvals(nonsig_subj_count,:) = {cluster_pvals}; %matrix of pvals corresponding to pixels in clusters (empty if none)
            subj_nonsig_cluster_sig_pixel_prop(nonsig_subj_count,:) = {cluster_sig_pixel_prop}; %proportion of cluster pvalues <= 0.05 (zero if no good clust)
            all_adj_pval_mats(e,:) = {adj_p};
            
            
            
        end 
        
        close all %close figure window 
        % just lots of if statements, no inner loops
        
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
    sig_elecs_by_subj(5,s) = {subj_sig_norm_plv_mats}; %norm plv mats for each sig elec
    sig_elecs_by_subj(6,s) = {subj_sig_adj_pvals}; %adjusted pvals from fdr correction for each sig elec
    sig_elecs_by_subj(7,s) = {subj_sig_clust_amps}; % sig elec amp freq index (corrected!)
    sig_elecs_by_subj(8,s) = {subj_sig_clust_phases}; % sig elec phase freq index
    sig_elecs_by_subj(9,s) = {subj_sig_clust_mats}; % sig clusters norm plv for each elec 
    sig_elecs_by_subj(10,s) = {subj_sig_clust_pvals}; % sig clusters pixels pvals for each elec
    sig_elecs_by_subj(11,s) = {subj_sig_cluster_sig_pixel_prop}; % prop of cluster pvalues <=0.05 for each elec (should all be >=0.95)


    nonsig_elecs_by_subj(1,s) = {nonsig_subj_count}; %num nonsig elecs within subj
    nonsig_elecs_by_subj(2,s) = {subj_nonsig_elecs}; %elec id of nonsig elecs within subj
    nonsig_elecs_by_subj(3,s) = {nonsig_mean_raw_plv}; %all_elec_mean_raw_plv
    nonsig_elecs_by_subj(4,s) = {nonsig_mean_norm_plv}; %all_elec_mean_norm_plv
    nonsig_elecs_by_subj(5,s) = {subj_nonsig_norm_plv_mats}; %norm plv mats for each sig elec
    nonsig_elecs_by_subj(6,s) = {subj_nonsig_adj_pvals}; %adjusted pvals from fdr correction for each nonsig elec
    nonsig_elecs_by_subj(7,s) = {subj_nonsig_clust_amps}; %amp freqs of largest cluster for nonsig elecs (zero if no clusters)
    nonsig_elecs_by_subj(8,s) = {subj_nonsig_clust_phases}; %phase freqs of largest cluster for nonsig elecs (zero if no clusters)
    nonsig_elecs_by_subj(9,s) = {subj_nonsig_clust_mats}; % clusters norm plv mats for each non elec (zero if no clusters)
    nonsig_elecs_by_subj(10,s) = {subj_nonsig_clust_pvals}; % clusters pixels pvals for each non elec (zero if no clusters)
    nonsig_elecs_by_subj(11,s) = {subj_nonsig_cluster_sig_pixel_prop}; % prop of cluster pvalues <=0.05 for each elec (should all be <0.95 or 0)
    
    %data to update all subj struct  
    %add correction info 
    data.fdr_adj_pvals = all_adj_pval_mats; %fdr corrected pval mats
    
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


end



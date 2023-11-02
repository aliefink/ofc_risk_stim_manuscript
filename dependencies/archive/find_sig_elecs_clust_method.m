function [sig_elecs_by_subj,nonsig_elecs_by_subj,total_count,nonsig_count,cfc_all_subj] = find_sig_elecs_clust_method(cfc_all_subj,subjs,pval_type)
%pval_type either 'prop' or 't_stat'

sig_elecs_by_subj = {};%1=num_sig_elecs; 2 sig elec ids; 3 num sig couples; 4 norm plv mats; 5 mean_plv plv; 6 adjusted pvals
nonsig_elecs_by_subj = {};%1=num_sig_elecs; 2 sig elec ids; 3 num sig couples; 4 norm plv mats; 5 mean_plv plv; 6 adjusted pvals


total_count = 0;
nonsig_count = 0;

for s=1:length(subjs)   
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    

    subj_count=0; %num significant electrodes within subj
    subj_sig_elecs = []; %elec ids of subj sig electrodes
    subj_prop_sig = []; 
    subj_sig_couples = []; %number of significant couplings for each significant electrode
    subj_sig_norm_plv = {}; %zscored plv matrix
    subj_sig_mean_norm_plv = []; %mean_plv of normalized plv for that electrode
    adj_plv_pvals = {}; %adjusted pvals from fdr correction
    
    
    nonsig_subj_count=0; %num significant electrodes within subj
    subj_nonsig_elecs = []; %elec ids of subj sig electrodes
    subj_prop_nonsig = [];
    subj_nonsig_couples = []; %number of significant couplings for each significant electrode
    subj_nonsig_norm_plv = {}; %zscored plv matrix
    subj_nonsig_mean_norm_plv = []; %mean_plv of normalized plv for that electrode
    adj_nonsig_pvals = {}; %adjusted pvals from fdr correction
    
    for e=1:data.n_elecs
        elec = data.ofc_elecs(e);
        mean_plv = data.mean_norm_plvs{e};
        plv = data.norm_real_plvs{e};
        if strcmp(pval_type,'prop') 
            norm_plv_pvals = data.plv_pvals_all{e}; %pvalues calculated from perm for each pixel (#perm > pixel/# perms total)
        elseif strcmp(pval_type, 't_stat')
            norm_plv_pvals = data.plv_tstat_all{e}; %pvalues calculated from perm for each pixel (t test)
        end 
% Benjamani & Hochberg False Discovery Rate (FDR) Correction
            % fdr_bh, 2015, David Groppe
        % must determine if pvalues are significant on *group-level*
         % within electrode - control for multiple comparisons
        %input pval matrix, false discovery rate q (if not 0.05), method =
        %pdep or dep, report = 'yes' or 'no'
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(norm_plv_pvals,0.05,'dep');
        n_sig_pixels = sum(h,'all'); %h = binary vector, h=1 is significant
        prop_sig_pixels = n_sig_pixels/760;
        
        %get contours of zscored matrix
        [cm,h] = contourf(plv);
        [contourTable, contourArray] = getContourLineCoordinates(cm);
%         if max(contourTable.Level) >= 4
%             max_z = 4; %if patient goes up to 5 std above mean - all levels above that are counted in significant contour 
%             
%         else 
%             max_z = 2;
%         end 
        if max(contourTable.Level) <= 2 %testing max levels 
            max_z = 2; 
%             max_level = contourTable.Level >= max_z; %corresponds to p<0.05 if =2
        elseif (max(contourTable.Level) > 2 && max(contourTable.Level)<=4)
            max_z = (max(contourTable.Level)-1); %setting max to max cluster + second largest cluster
%             max_level = contourTable.Level >= max_z; %extracts all values of 2 highest level clusters 
        elseif (max(contourTable.Level) > 4 && max(contourTable.Level)<=6)
            max_z = 4; %if patient goes up to 5 std above mean - all levels above that are counted in significant contour 
%             max_level = contourTable.Level>=max_z;
        elseif max(contourTable.Level) > 6
            max_z = 6; %if patient goes up to 5 std above mean - all levels above that are counted in significant contour 
            
        end 
        max_level = contourTable.Level>=max_z;

        if sum(max_level) == 0
            prop_sig_clust = 0;
        else 
            max_levelX = contourTable.X(max_level);
            max_levelY = contourTable.Y(max_level); %amplitudes of most significant cluster
            amp_f_array = 5:5:200;
            phase_f_array = 2:20;
            sig_phase_idx = unique(round(max_levelX));
            sig_amp_idx = unique(round(max_levelY));
            sig_phases = phase_f_array(sig_phase_idx');
            sig_amps = amp_f_array(sig_amp_idx');
            cluster_mat = plv(sig_amp_idx,sig_phase_idx);
            numel(cluster_mat); %testing cluster size!
            cluster_pvals = adj_p(sig_amp_idx,sig_phase_idx);
            prop_sig_clust = sum(cluster_pvals<=0.05,'all')/numel(cluster_pvals);
        end
        

        if prop_sig_clust >=0.95

            subj_count=subj_count+1; %add one to num sig elecs within subj
            total_count = total_count+1; %add one to total num sig elecs
            subj_sig_elecs(subj_count,:) = elec; %add elec id to subj tracker
            subj_prop_sig(subj_count,:) = prop_sig_pixels;
            subj_sig_couples(subj_count,:) = n_sig_pixels; %add num sig couples to subj tracker
            subj_sig_norm_plv(subj_count,:) = {plv};
            subj_sig_mean_norm_plv(subj_count,:) = mean_plv;
            adj_plv_pvals(subj_count,:) = {adj_p};
            
            
        else
            nonsig_subj_count = nonsig_subj_count+1;
            nonsig_count = nonsig_count+1;
            subj_nonsig_elecs(nonsig_subj_count,:) = elec;
            subj_prop_nonsig(nonsig_subj_count,:) = prop_sig_pixels;
            subj_nonsig_couples(nonsig_subj_count,:) = n_sig_pixels;
            subj_nonsig_norm_plv(nonsig_subj_count,:) = {plv};
            subj_nonsig_mean_norm_plv(nonsig_subj_count,:) = mean_plv;
            adj_nonsig_pvals(nonsig_subj_count,:) = {adj_p};
        end 
        
        
    end 

    sig_elecs_by_subj(1,s) = {subj_count}; 
    sig_elecs_by_subj(2,s) = {subj_sig_elecs};
    sig_elecs_by_subj(3,s) = {subj_sig_couples};
    sig_elecs_by_subj(4,s) = {subj_sig_norm_plv};
    sig_elecs_by_subj(5,s) = {subj_sig_mean_norm_plv};
    sig_elecs_by_subj(6,s) = {adj_plv_pvals};
    
    nonsig_elecs_by_subj(1,s) = {nonsig_subj_count}; 
    nonsig_elecs_by_subj(2,s) = {subj_nonsig_elecs};
    nonsig_elecs_by_subj(3,s) = {subj_nonsig_couples};
    nonsig_elecs_by_subj(4,s) = {subj_nonsig_norm_plv};
    nonsig_elecs_by_subj(5,s) = {subj_nonsig_mean_norm_plv};
    nonsig_elecs_by_subj(6,s) = {adj_nonsig_pvals};
    
    data.sig_elecs = subj_sig_elecs;
    mean_sig_plv = mean(subj_sig_mean_norm_plv);
    data.mean_plv_sig_plv = mean_sig_plv;
    %add non significant data to cfc_all_subj;
    data.nonsig_elecs = subj_nonsig_elecs;
    mean_nonsig_plv = mean(subj_nonsig_mean_norm_plv);
    data.mean_plv_nonsig_plv = mean_nonsig_plv;
    cfc_all_subj.(subj_id) = data;

    

    

end


end


%noncluster method
% n_sig_pixels = sum(h,'all'); %h = binary vector, h=1 is significant
% prop_sig_pixels = n_sig_pixels/760;
%         if max(contourTable.Level) >= 5
%             max_z = 5; %if patient goes up to 5 std above mean - all levels above that are counted in significant contour 
%             
%         else 
%             max_z = 2;

%% clustering archive 
        %         clust_idx = idx(corepts);
%         clust_ids = unique(clust_idx);
%         num_clust = length(clust_ids);
%         clust_p = [];
%         for c=1:range(num_clust)
%             c_num = clust_ids(c);
%         end 
%         
%         if n_sig_pixels ~= 0
%         mean_p = mean(adj_p,'all');
%         if clust_alpha <= 0.05
%         if prop_sig_pixels >= .95 % 10% are significant use this for prop!


         %phases of most significant cluster
        

% 
        % clustering attempt lol
%         [idx,V,D] = spectralcluster(adj_p,40);
%         max_clust = mode(idx); %find index of largest cluster 
%         max_clust_data = adj_p(find(idx==max_clust));
%         max_clust_data(max_clust_data>1) = 1;
%         clust_alpha = (sum(max_clust_data(find(idx==max_clust))<0.05)/length(max_clust_data(find(idx==max_clust))));
%         minpts = 4;%recommended to set equal to num of dimensions of input data+1 (phase freqs+1)
%         [idx,corepts] = dbscan(adj_p,1,minpts); %test for significant clusters in sample 
%         number_of_clusters = sum(unique(idx)>0);
%         max_clust = corepts(find(idx==1));
%         core_data = adj_p(find(idx==1), 1);
%         core_data(core_data>1) = 1; 
%         core_idx = idx(corepts);
%         clusters = splitapply(@(x) {x}, core_data, core_idx);
%         max_clust_data = core_data(core_idx);
%         prop_sig_clust = sum(core_data<0.05,'all')/numel(core_data);
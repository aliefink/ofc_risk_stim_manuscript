function [cfc_all_subj,all_elec_summary_stats]...
    = get_summary_stats(cfc_all_subj, subjs)


%for each subject I need - 
    %max phase freq - per electrode & subj total
    %max amp freq - per electrode & subj total
    %max sig plv - per electrode & subj total
    

phase_f_array = 2:20;
amp_f_array = 5:5:200;
    
for s=1:length(subjs) %iterate through all subjects   
    subj_id = char(subjs(s));
    data = cfc_all_subj.(subj_id);
    
    subj_sig_pixels = {};
    subj_max_phase_by_elec = {};
    subj_max_amp_by_elec = {};
    
    for e=1:data.n_elecs
        
        elec = data.ofc_elecs(e);
        plv_z = data.norm_plvs{e};
        adj_p = data.fdr_adj_pvals{e,1};
        [row,col] = find(adj_p < 0.05); %find indices of sig plvs (corrected pval < 0.05)
        sig_plv = plv_z(row,col); %extract norm plv values for sig pixels
        adj_plv = zeros(size(plv_z)); %make empty adj plv matrix
        adj_plv(row,col) = sig_plv; %include plvs for only sig pixels
        subj_sig_pixels(e) = {adj_plv};
        phase_means = mean(adj_plv,1);
        max_phase_plv = max(phase_means);
        phase_max = phase_f_array(find(phase_means==max_phase_plv));
        subj_max_phase_by_elec(e) = {phase_max};
        
        amp_means = mean(adj_plv,2)'; %mean along rows - amp means
        max_amp_plv = max(amp_means);
        amp_max = amp_f_array(find(amp_means==max_amp_plv));
        subj_max_amp_by_elec(e) = {amp_max};
    
    
    end 
    
    data.sig_pixels = subj_sig_pixels;
    data.max_phases = subj_max_phase_by_elec;
    data.max_amps = subj_max_amp_by_elec;
    
    % update cfc_all_subj subj field updated data structure
    cfc_all_subj.(subj_id) = data;
    
    

end
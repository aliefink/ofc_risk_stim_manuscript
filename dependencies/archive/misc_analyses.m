%% RANDOM PLOTS

% %get plot data
% %         
% % [all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
% %     get_plv_freq_block_method(cfc_all_subj, subjs);
% % 
% % freq_specific_plvs = {all_elec_plv_info.delta_theta_gamma_means, all_elec_plv_info.delta_theta_hga_means,...
% %     all_elec_plv_info.alpha_gamma_means, all_elec_plv_info.alpha_hga_means,...
% %     all_elec_plv_info.beta_gamma_means,all_elec_plv_info.beta_hga_means};
% % hga_letters = "HG";
% 
% 
% plv_by_freq = figure('Name','Mean PLV by Freq');
% plv_by_freq.Position = [0,0,750, 600];% plv_by_freq.Position(3:4) = [1000 1000];
% 
% x_idxs = [1 2 3];
% bar_y_vals = [mean(all_elec_plv_info.delta_theta_gamma_means), mean(all_elec_plv_info.delta_theta_hga_means);...
%     mean(all_elec_plv_info.alpha_gamma_means), mean(all_elec_plv_info.alpha_hga_means);...
%     mean(all_elec_plv_info.beta_gamma_means),mean(all_elec_plv_info.beta_hga_means)];
% b = bar(x_idxs,bar_y_vals,'grouped','FaceColor','flat','EdgeColor', 'w');
% set(gca,'XTickLabels', {'\fontname{Lucida Grande}Delta-Theta \bf\delta\theta\rm',...
%     '\fontname{Lucida Grande}Alpha \bf\alpha\rm','\fontname{Lucida Grande}Beta \bf\beta\rm'},...
%     'FontSize',22,'FontName','Lucida Grande')%...
% 
% title({'\rm\fontname{Lucida Grande}Low Frequency Phase Coupling to Gamma \fontsize{34}(\bf\gamma\rm)', '\fontsize{30}& High Gamma (HG) Amplitudes'},...
%     'VerticalAlignment','bottom','FontSize',30)
% xlabel('\fontname{Lucida Grande}Frequency for Phase (Hz)','FontSize',26)
% ylabel('\fontname{Lucida Grande}Phase-Amplitude Coupling (PLV)','FontSize',26)
% hold on 
% %manually set color data for each bar
% b(1,1).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
% b(1,1).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
% b(1,1).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;[252,174,145]./255;0.9961    0.8980    0.8510
% 
% b(1,2).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
% b(1,2).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
% b(1,2).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;0.9957,0.8784,0.8235
% 
% %add amp frequency label to top of bars
% xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
% ytips1 = b(1,1).YEndPoints; %1x3 Y end points of gamma vals
% labels1 = {'\fontsize{20}\fontname{Lucida Grande}\bf\delta\theta\cdot\fontsize{22}\gamma',...
%     '\fontsize{20}\fontname{Lucida Grande}\bf\alpha\cdot\fontsize{22}\gamma',...
%     '\fontsize{20}\fontname{Lucida Grande}\bf\beta\cdot\fontsize{22}\gamma'};%should be same length as xtips (1x3)
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% xtips2 = b(1,2).XEndPoints; %1x3 X end points of hga vals
% ytips2 = b(1,2).YEndPoints; %1x3 Y end points of hga vals
% labels2 = {join(['\fontname{Lucida Grande}\fontsize{20}\bf\delta\theta\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2),...
%     join(['\fontname{Lucida Grande}\fontsize{20}\bf\alpha\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2),...
%     join(['\fontname{Lucida Grande}\fontsize{20}\bf\beta\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2)};%should be same length as xtips (1x3)
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% 
% % title({'\rm\fontname{Lucida Grande}Low Frequency Phase Coupling with High Frequency Amplitudes:'; 'Gamma (\bf\gamma\rm) Oscillations & High Frequency Activity (\itHFA\rm)'},...
% %     'VerticalAlignment','bottom','FontSize',24)
% % 
% % set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% % % set(gcf,'PaperPosition', [0 0 20 10]); %set position of fig to be consistent with pdf sizing
% % print(gcf,[fig_path 'barplot.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page
% 
% % alternative saving options:
% % exportgraphics(gcf,[fig_path 'barplot_plv_by_freq.pdf'],"ContentType","vector")
% % saveas(plv_by_freq,[fig_path 'barplot_test.pdf'])
% 
% % set(gca,'XTickLabels', {'\fontsize{20}\fontname{Lucida Grande}\bf\delta\theta\cdot\fontsize{22}\gamma',...
% %     '\fontsize{20}\fontname{Lucida Grande}\bf\alpha\cdot\fontsize{22}\gamma',...
% %     '\fontsize{20}\fontname{Lucida Grande}\bf\beta\cdot\fontsize{22}\gamma'},...
% %     'FontSize',22,'FontName','Lucida Grande')%...
% %% BROADBAND GAMMA ONLY 
% 
% plv_by_freq = figure('Name','Mean PLV by Freq');
% plv_by_freq.Position = [0,0,700, 600];% plv_by_freq.Position(3:4) = [1000 1000];
% 
% x_idxs = [1 2 3];
% bar_y_vals = [mean(all_elec_plv_info.delta_theta_broadgamma_means);...
%     mean(all_elec_plv_info.alpha_broadgamma_means);...
%     mean(all_elec_plv_info.beta_broadgamma_means)];
% 
% b = bar(x_idxs,bar_y_vals,'FaceColor','flat','EdgeColor', 'w');
% 
% set(gca,'XTickLabels', {'\fontname{Lucida Grande}Delta-Theta',...
%     '\fontname{Lucida Grande}Alpha',...
%     '\fontname{Lucida Grande}Beta'},...
%     'FontSize',22,'FontName','Lucida Grande')
% 
% title({'\rm\fontname{Lucida Grande}Coupling Between Low Frequency Phase and'; 'Broadband Gamma Amplitude'},...
%     'VerticalAlignment','bottom','FontSize',24)
% 
% 
% xlabel('\fontname{Lucida Grande}Frequency for Phase (Hz)','FontSize',24)
% ylabel('\fontname{Lucida Grande}Phase-Amplitude Coupling (PLV)','FontSize',24)
% 
% hold on 
% 
% 
% %manually set color data for each bar
% b(1,1).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
% b(1,1).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
% b(1,1).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;
% 
% hold on 
% %add sem bars to amp frequency label to top of bars
% %calculate sem for bar
% b1_data_pts = all_elec_plv_info.delta_theta_broadgamma_means;
% b1_sem =  std(all_elec_plv_info.delta_theta_broadgamma_means)/sqrt(15);
% b1_upper = mean(b1_data_pts)+b1_sem;
% b1_lower = mean(b1_data_pts)-b1_sem;
% 
% b2_data_pts = all_elec_plv_info.alpha_broadgamma_means;
% b2_sem =  std(all_elec_plv_info.alpha_broadgamma_means)/sqrt(15);
% b2_upper = mean(b2_data_pts)+b2_sem;
% b2_lower = mean(b2_data_pts)-b2_sem;
% 
% b3_data_pts = all_elec_plv_info.beta_broadgamma_means;
% b3_sem =  std(all_elec_plv_info.beta_broadgamma_means)/sqrt(15);
% b3_upper = mean(b3_data_pts)+b3_sem;
% b3_lower = mean(b3_data_pts)-b3_sem;
% 
% 
% % err_data = [mean(all_elec_plv_info.delta_theta_broadgamma_means),...
% %     mean(all_elec_plv_info.alpha_broadgamma_means),...
% %     mean(all_elec_plv_info.beta_broadgamma_means)];
% err_data = bar_y_vals;
% err_sem = [b1_sem;b2_sem;b3_sem];
% % [ngroups, nbars] = size(bar_y_vals);
% % groupwidth = min(0.8, nbars/(nbars + 1.5));
% % x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% 
% %get x coords of bars
% xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
% errorbar(xtips1',err_data,err_sem,'k','linestyle','none','LineWidth', 1.5);%'Color', [.25 .25 .25]);
% % err_lower = [b1_lower,b2_lower,b3_lower];
% % err_upper = [b1_upper,b2_upper,b3_upper];
% 
% 
% % xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
% % ytips1 = b(1,1).YEndPoints; %1x3 Y end points of gamma vals
% 
% 
% % labels1 = {'\fontname{Lucida Grande}\bf\delta\theta\cdot\gamma',...
% %     '\fontname{Lucida Grande}\bf\alpha\cdot\gamma',...
% %     '\fontname{Lucida Grande}\bf\beta\cdot\gamma'};%should be same length as xtips (1x3)
% % text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
% %     'VerticalAlignment','bottom','FontSize',18)
% % set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% % % set(gcf,'PaperPosition', [0 0 20 10]); %set position of fig to be consistent with pdf sizing
% % print(gcf,[fig_path 'barplot_broadgamma_sem.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page



% %% Plot percent sig elecs within subj x risk pref 
% risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% % percent_sig_elecs = [0 34.69387755 85.71428571 9.090909091 70 42.85714286 33.33333333 0 33.33333333 50 33.33333333 0 0 0 0]; %prop
% percent_sig_elecs = [0 30.6122449 85.71428571 9.090909091 70 14.28571429 0 0 33.33333333 50 33.33333333 0 0 0 0]; %tstat
% % percent_sig_elecs = [20 44.89795918 57.14285714 72.72727273 50 42.857142860 0 33.33333333 50 33.33333333 0 0 16.66666667 0]; - gamma distributioncut offs
% 
% 
% lm = fitlm(percent_sig_elecs,risk_prefs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(risk_prefs.',percent_sig_elecs.','type','pearson');
% 
% scatter(percent_sig_elecs,risk_prefs)
% xlabel('Percent of Electrodes with Significant Coupling by Subject')
% ylabel('Risk Preference')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))



% %% Average PLVs for NONsignificant electrodes
% 
% risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% average_nonsig_plvs = [];
% 
% 
% for s=1:length(subjs)
%     subj_id = char(subjs(s));
%     data = nonsig_elecs_by_subj{4,s}; %mean plv values for significant elecs
%     mean_nonsig = mean(data);
%     average_nonsig_plvs(1,s) = mean_nonsig;
%     
% end 
% 
% 
% average_nonsig_plvs(isnan(average_nonsig_plvs)) = 0;
% 
% lm = fitlm(average_nonsig_plvs,risk_prefs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(risk_prefs.',average_nonsig_plvs.','type','pearson');
% 
% scatter(average_nonsig_plvs,risk_prefs)
% xlabel('Average PLV of Non Significant Electrodes by Subject')
% ylabel('Risk Preference')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))


% %% Plot average Significant cluster within subj x risk pref
% 
% risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% clust_cfc_bysubj = [];
% 
% for s=1:length(subjs)
%     subj_id = char(subjs(s));
%     data = cfc_all_subj.(subj_id);
%     clust_cfc_means = [];
%     for e=1:data.n_elecs
%         plv = data.norm_real_plvs{1,e};
%         clust_cfc_plv = plv(sig_amp_idx',sig_phase_idx');
%         clust_cfc_means(1,e) = mean(clust_cfc_plv,'all');
%     end 
%     clust_cfc_bysubj(s) = mean(clust_cfc_means);
% end 
% 
% 
% 
% lm = fitglm(clust_cfc_bysubj,risk_prefs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(risk_prefs.',clust_cfc_bysubj.','type','pearson');
% 
% scatter(clust_cfc_bysubj,risk_prefs)
% xlabel('Average Cluster CFC PLV of Significant Electrodes by Subject')
% ylabel('Risk Preference')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))
% 


% 
% %% Plot average BETA-DELTA/THETA of Significant electrodes within subj x risk pref
% 
% risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% beta_cfc_bysubj = [];
% 
% for s=1:length(subjs)
%     subj_id = char(subjs(s));
%     data = cfc_all_subj.(subj_id);
%     beta_theta_cfc_means = [];
%     for e=1:data.n_elecs
%         plv = data.norm_real_plvs{1,e};
%         beta_cfc_plv = plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
%         beta_theta_cfc_means(1,e) = mean(beta_cfc_plv,'all');
%     end 
%     beta_cfc_bysubj(s) = mean(beta_theta_cfc_means);
% end 
% 
% 
% 
% lm = fitglm(beta_cfc_bysubj,risk_prefs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(risk_prefs.',beta_cfc_bysubj.','type','pearson');
% 
% scatter(beta_cfc_bysubj,risk_prefs)
% xlabel('Average beta CFC PLV of Significant Electrodes by Subject')
% ylabel('Risk Preference')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))
% 

% 
% %% Plot average GAMMA-DELTA/THETA of Significant electrodes within subj x risk pref
% 
% risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% gamma_cfc_bysubj = [];
% 
% for s=1:length(subjs)
%     subj_id = char(subjs(s));
%     data = cfc_all_subj.(subj_id);
%     gamma_theta_cfc_means = [];
%     for e=1:data.n_elecs
%         plv = data.norm_real_plvs{1,e};
%         gamma_cfc_plv = plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
%         gamma_theta_cfc_means(1,e) = mean(gamma_cfc_plv,'all');
%     end 
%     gamma_cfc_bysubj(s) = mean(gamma_theta_cfc_means);
% end 
% 
% 
% 
% lm = fitglm(gamma_cfc_bysubj,risk_prefs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(risk_prefs.',gamma_cfc_bysubj.','type','pearson');
% 
% scatter(gamma_cfc_bysubj,risk_prefs)
% xlabel('Average gamma CFC PLV of Significant Electrodes by Subject')
% ylabel('Risk Preference')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))
% 


% %% Plot average HGA-DELTA/THETA of Significant electrodes within subj x risk pref
% risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% hga_cfc_bysubj = [];
% 
% 
% for s=1:length(subjs)
%     subj_id = char(subjs(s));
%     data = cfc_all_subj.(subj_id);
%     hga_theta_cfc_means = [];
%     for e=1:data.n_elecs
%         plv = data.norm_real_plvs{1,e};
%         hga_cfc_plv = plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
%         hga_theta_cfc_means(1,e) = mean(hga_cfc_plv,'all');
%     end 
%     hga_cfc_bysubj(s) = mean(hga_theta_cfc_means);
% end 
% 
% 
% 
% lm = fitglm(hga_cfc_bysubj,risk_prefs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(risk_prefs.',hga_cfc_bysubj.','type','pearson');
% 
% scatter(hga_cfc_bysubj,risk_prefs)
% xlabel('Average HGA CFC PLV of Significant Electrodes by Subject')
% ylabel('Risk Preference')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))
% 

% 
% 
% 
% %% Plot average PLV of Significant electrodes within subj x risk pref
% 
% risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% average_sig_plvs = [];
% 
% 
% for s=1:length(subjs)
%     subj_id = char(subjs(s));
%     data = sig_elecs_by_subj{5,s}; %mean plv values for significant elecs
%     mean_sig = mean(data);
%     average_sig_plvs(1,s) = mean_sig;
%     
% end 
% 
% 
% average_sig_plvs(isnan(average_sig_plvs)) = 0;
% 
% lm = fitglm(average_sig_plvs,risk_prefs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(risk_prefs.',average_sig_plvs.','type','pearson');
% 
% scatter(average_sig_plvs,risk_prefs)
% xlabel('Average PLV of Significant Electrodes by Subject')
% ylabel('Risk Preference')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))
% 
% 
% 
% %% Plot average PLV of ALL electrodes within subj x risk pref 
% 
% risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
% average_subj_plvs = [];
% 
% for s=1:length(subjs)
%     subj_id = char(subjs(s));
%     data = cfc_all_subj.(subj_id);
%     mean_plv_all_elec = mean(cell2mat(data.mean_norm_plvs));
%     average_subj_plvs(1,s) = mean_plv_all_elec;
%     
% 
% end 
% 
% 
% 
% lm = fitglm(average_subj_plvs,risk_prefs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(risk_prefs.',average_subj_plvs.','type','pearson');
% 
% scatter(average_subj_plvs,risk_prefs)
% xlabel('Average PLV of All Electrodes by Subject')
% ylabel('Risk Preference')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))
% 
% 
% 
% %% Supplementary Figures 
% 
% % Histogram of permutation values 
% 
% 
% surr = squeeze(cfc_all_subj.s02.perm_plvs{1,1}(1,1,:));
% 
% plv_zscore = cfc_all_subj.s02.norm_real_plvs{1,1}(1,1);
% 
% perm_fig = figure('Name','Single Electrode Single Pixel Permutation Distribution')
% histogram(surr,20,'FaceColor','#D95319')
% hold on 
% xline(plv_zscore,'--k','Observed PLV zscore','LineWidth',1.5)
% ylabel('Frequency')
% xlabel('Surrogate Phase Locking Values')
% title('Example Significant Permutation Distribution (n=1000)')
% saveas(perm_fig,[fig_path 'permutation_example.pdf'])
% 
% % Another way to calculate pvalue - gives slightly different result...
% % F = @(x)(exp (-0.5*(x.^2))./sqrt (2*pi));
% % p = quad(F, plv_zscore, 100);
% 
% 
% %% Figures single elec plv 
% fig_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/figs/cfc/results_methods/';
% addpath(genpath(fig_path))
% 
% data = cfc_all_subj.s06;
% 
% elec_num = 11;
% plv_mat = data.norm_real_plvs{1,find(data.ofc_elecs==elec_num)};
% plv_single_elec = figure('Name','zscore PLV one subj one elec')
% contourf(plv_mat);
% amp_f_array = 5:5:200;
% phase_f_array = 2:20;
% % Set locations of ticks and their labels.
% set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
% 'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',14);
% title(['Mean Phase Amplitude Coupling'],'FontSize',16);
% xlabel('frequency for phase (Hz)','FontSize',14);
% ylabel('frequency for amplitude (Hz)','FontSize',14);
% h = colorbar;
% set(get(h, 'ylabel'), 'string', 'PLV Zscore','FontSize',14);
% saveas(plv_single_elec,[fig_path, num2str(elec_num), '_single_elec_plv.pdf'])
% 
% 
% 
% % plot mean amp ranges 
% % 
% % [sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, sig_elecs_by_subj, nonsig_elecs_by_subj, subjs, total_count,nonsig_count);
% 
% 
% sig_plv_group = figure('Name','mean zscore PLV significant elecs across subj')
% contourf(sig_mean_plv);
% amp_f_array = 5:5:200;
% phase_f_array = 2:20;
% % Set locations of ticks and their labels.
% set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
% 'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',14);
% title(['Mean Phase Amplitude Coupling'],'FontSize',16);
% xlabel('frequency for phase (Hz)','FontSize',14);
% ylabel('frequency for amplitude (Hz)','FontSize',14);
% h = colorbar;
% set(get(h, 'ylabel'), 'string', 'PLV Zscore','FontSize',14);
% saveas(sig_plv_group,[fig_path 'group_mean_plv_sig_elecs.pdf'])
% 
% 




%% RANDOM ANALYSES

% %% extract contours from contourf plot results 
% 
% [cm, h] = contourf(sig_mean_plv);
% % contourTable = getContourLineCoordinates(cm);
% [contourTable, contourArray] = getContourLineCoordinates(cm);
% 
% % get X,Y coords of most significant contours (> 3 std)
% % max_level = contourTable.Level; == h.LevelList(7) | contourTable.Level == h.LevelList(8);
% max_level = contourTable.Level > 3;
% max_levelX = contourTable.X(max_level); %phases of most significant cluster
% max_levelY = contourTable.Y(max_level); %amplitudes of most significant cluster
% 
% amp_f_array = 5:5:200;
% phase_f_array = 2:20;
% 
% sig_phase_idx = unique(round(max_levelX));
% sig_amp_idx = unique(round(max_levelY));
% sig_phases = phase_f_array(sig_phase_idx');
% sig_amps = amp_f_array(sig_amp_idx');
% 
% % Getting Mean PLV by freq band
% 
% %define frequency bands
% delta_theta = [2 8];
% alpha = [9 12];
% beta = [13 29];
% gamma = [30 60];
% hga = [61 200];
% 
% %find freqs by significant phases
% delta_theta_phases = sig_phases(sig_phases>= delta_theta(1)&sig_phases<= delta_theta(2));
% alpha_phases = sig_phases(sig_phases>= alpha(1)&sig_phases<= alpha(2));
% 
% %find freqs by significant amps
% alpha_amps = sig_amps(sig_amps>= alpha(1)&sig_amps<= alpha(2));
% beta_amps = sig_amps(sig_amps>= beta(1)&sig_amps<= beta(2));
% gamma_amps = sig_amps(sig_amps>= gamma(1)&sig_amps<= gamma(2));
% hga_amps = sig_amps(sig_amps>= hga(1)&sig_amps<= hga(2));
% 
% 
% % statistical test between beta amp x delta-theta phase coupling
% beta_delta_theta_amp_idx = find(ismember(amp_f_array,beta_amps));
% beta_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_beta_delta_theta_mean_cfc = sig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
% nonsig_beta_delta_theta_mean_cfc = nonsig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
% 
% [beta_delta_theta_clusters, beta_delta_theta_p_values, beta_delta_theta_t_sums, beta_delta_theta_permutation_distribution ] = permutest(sig_beta_delta_theta_mean_cfc,nonsig_beta_delta_theta_mean_cfc);
% disp(beta_delta_theta_p_values)
% 
% mean_beta_delta_theta_cfc = mean(sig_beta_delta_theta_mean_cfc,'all');
% 
% 
% % statistical test between gamma amp x delta-theta phase coupling
% gamma_delta_theta_amp_idx = find(ismember(amp_f_array,gamma_amps));
% gamma_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_gamma_delta_theta_mean_cfc = sig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
% nonsig_gamma_delta_theta_mean_cfc = nonsig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
% 
% [gamma_delta_theta_clusters, gamma_delta_theta_p_values, gamma_delta_theta_t_sums, gamma_delta_theta_permutation_distribution ] = permutest(sig_gamma_delta_theta_mean_cfc,nonsig_gamma_delta_theta_mean_cfc);
% disp(gamma_delta_theta_p_values)
% mean_gamma_delta_theta_cfc = mean(sig_gamma_delta_theta_mean_cfc,'all');
% 
% 
% 
% % statistical test between hga amp x delta-theta phase coupling
% hga_delta_theta_amp_idx = find(ismember(amp_f_array,hga_amps));
% hga_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_hga_delta_theta_mean_cfc = sig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
% nonsig_hga_delta_theta_mean_cfc = nonsig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
% 
% [hga_delta_theta_clusters, hga_delta_theta_p_values, hga_delta_theta_t_sums, hga_delta_theta_permutation_distribution ] = permutest(sig_hga_delta_theta_mean_cfc,nonsig_hga_delta_theta_mean_cfc);
% disp(hga_delta_theta_p_values)
% mean_hga_delta_theta_cfc = mean(sig_hga_delta_theta_mean_cfc,'all');
% 
% 


%%%% extract the frequencies of amp/phase that are in the most significant
%%%% contours

cm = contourf(sig_mean_plv);
close all
contourTable = getContourLineCoordinates(cm);

% get indices from most sig clusters
max_z_all = max(contourTable.Level); %max level mean plv all sig elecs
sig_thresh = 2; %setting threshold for 2 std above mean - level > 2 *may need to raise (p < 0.05)
max_level = contourTable.Level > sig_thresh; %corresponds to p<0.001
max_levelX = contourTable.X(max_level); %extract x index of max contours (phase freq)
max_levelY = contourTable.Y(max_level); %extract y index of max contours (amp freq)

% get freq info from clust indices (must round contour index to whole number - actual data indices can't be exact contour curve)
sig_amp_idx = unique(round(max_levelY)); %indexes of contour amp freqs (rounded)
sig_phase_idx = unique(round(max_levelX)); %indexes of contour phase freqs (rounded)
amp_f_array = 5:5:200; % values of amp freqs
phase_f_array = 2:20; % values of phase freqs
sig_amps = amp_f_array(sig_amp_idx'); %amp freqs corresponding to amp index
sig_phases = phase_f_array(sig_phase_idx'); %phase freqs corresponding to phase index


[all_elec_plv_info,sig_elec_plv_info] = get_plv_freq_info_clust_method(cfc_all_subj, subjs, sig_amp_idx,sig_phase_idx);

% %define frequency bands
% delta_theta = [2 8];
% beta = [13 30];
% gamma = [30 60];
% hga = [60 200];
% 
% %find freqs by significant phases
% delta_theta_phases = sig_phases(sig_phases>= delta_theta(1)&sig_phases<= delta_theta(2));
% 
% %find freqs by significant amps
% beta_amps = sig_amps(sig_amps>= beta(1)&sig_amps< beta(2));
% gamma_amps = sig_amps(sig_amps>= gamma(1)&sig_amps< gamma(2));
% hga_amps = sig_amps(sig_amps>= hga(1)&sig_amps<= hga(2));
% 
% 
% 
% 
% 
% beta_delta_theta_amp_idx = find(ismember(amp_f_array,beta_amps));
% beta_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_beta_delta_theta_mean_cfc = sig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
% 
% gamma_delta_theta_amp_idx = find(ismember(amp_f_array,gamma_amps));
% gamma_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_gamma_delta_theta_mean_cfc = sig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
% 
% hga_delta_theta_amp_idx = find(ismember(amp_f_array,hga_amps));
% hga_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_hga_delta_theta_mean_cfc = sig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
% 

% 
% %% 
% 
% 
% 
% % statistical test between beta amp x delta-theta phase coupling
% beta_delta_theta_amp_idx = find(ismember(amp_f_array,beta_amps));
% beta_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_beta_delta_theta_mean_cfc = sig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
% nonsig_beta_delta_theta_mean_cfc = nonsig_mean_plv(beta_delta_theta_amp_idx,beta_delta_theta_phase_idx);
% 
% [beta_delta_theta_clusters, beta_delta_theta_p_values, beta_delta_theta_t_sums, beta_delta_theta_permutation_distribution ] = permutest(sig_beta_delta_theta_mean_cfc,nonsig_beta_delta_theta_mean_cfc);
% disp(beta_delta_theta_p_values)
% 
% mean_beta_delta_theta_cfc = mean(sig_beta_delta_theta_mean_cfc,'all');
% 
% 
% % statistical test between gamma amp x delta-theta phase coupling
% gamma_delta_theta_amp_idx = find(ismember(amp_f_array,gamma_amps));
% gamma_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_gamma_delta_theta_mean_cfc = sig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
% nonsig_gamma_delta_theta_mean_cfc = nonsig_mean_plv(gamma_delta_theta_amp_idx,gamma_delta_theta_phase_idx);
% 
% [gamma_delta_theta_clusters, gamma_delta_theta_p_values, gamma_delta_theta_t_sums, gamma_delta_theta_permutation_distribution ] = permutest(sig_gamma_delta_theta_mean_cfc,nonsig_gamma_delta_theta_mean_cfc);
% disp(gamma_delta_theta_p_values)
% mean_gamma_delta_theta_cfc = mean(sig_gamma_delta_theta_mean_cfc,'all');
% 
% 
% 
% % statistical test between hga amp x delta-theta phase coupling
% hga_delta_theta_amp_idx = find(ismember(amp_f_array,hga_amps));
% hga_delta_theta_phase_idx = find(ismember(phase_f_array,delta_theta_phases));
% sig_hga_delta_theta_mean_cfc = sig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
% nonsig_hga_delta_theta_mean_cfc = nonsig_mean_plv(hga_delta_theta_amp_idx,hga_delta_theta_phase_idx);
% 
% [hga_delta_theta_clusters, hga_delta_theta_p_values, hga_delta_theta_t_sums, hga_delta_theta_permutation_distribution ] = permutest(sig_hga_delta_theta_mean_cfc,nonsig_hga_delta_theta_mean_cfc);
% disp(hga_delta_theta_p_values)
% mean_hga_delta_theta_cfc = mean(sig_hga_delta_theta_mean_cfc,'all');

% 
% %% Get mean plv mat contours
% %check bwlabeln and bwconncomp for cluster extraction
% 
% %extract contours from mean zscore sig elecs plot
% cm = contourf(sig_mean_plv);
% close all
% contourTable = getContourLineCoordinates(cm);
% 
% % get indices from most sig clusters
% max_z_all = max(contourTable.Level); %max level mean plv all sig elecs
% sig_thresh = 2; %setting threshold for 2 std above mean - level > 2 *may need to raise (p < 0.05)
% max_level = contourTable.Level > sig_thresh; %corresponds to p<0.001
% max_levelX = contourTable.X(max_level); %extract x index of max contours (phase freq)
% max_levelY = contourTable.Y(max_level); %extract y index of max contours (amp freq)
% 
% % get freq info from clust indices (must round contour index to whole number - actual data indices can't be exact contour curve)
% sig_amp_idx = unique(round(max_levelY)); %indexes of contour amp freqs (rounded)
% sig_phase_idx = unique(round(max_levelX)); %indexes of contour phase freqs (rounded)
% amp_f_array = 5:5:200; % values of amp freqs
% phase_f_array = 2:20; % values of phase freqs
% sig_amps = amp_f_array(sig_amp_idx'); %amp freqs corresponding to amp index
% sig_phases = phase_f_array(sig_phase_idx'); %phase freqs corresponding to phase index

% %% Permutation Testing 

%permutation testing between mean plv matrices
% [mean_clusters, mean_p_values, mean_t_sums, mean_permutation_distribution ] = permutest(sig_mean_plv,nonsig_mean_plv);
% disp(mean_p_values);
% one_samp_plv_null = zeros(size(sig_mean_plv));
% [mean_clusters, mean_p_values, mean_t_sums, mean_permutation_distribution ] = permutest(sig_mean_plv,one_samp_plv_null);
% disp(mean_p_values);
% plv = cfc_all_subj.s02.real_plvs{1,1};
% perm_plv = mean(cfc_all_subj.s02.perm_plvs{1,1},3);
% [mean_clusters, mean_p_values, mean_t_sums, mean_permutation_distribution ] = permutest(plv,perm_plv);

%function [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(data,n_perm,tail,alpha_level,mu,reports,seed_state)

%permutation testing by specific frequency band

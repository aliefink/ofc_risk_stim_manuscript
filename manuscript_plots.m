%% data loading
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

subjs = {'s02','s04','s06','s07','s08','s09','s10','s12','s13','s14','s15','s16','s17','s18','s19'};

%% Data Loading


%cfc_all_subj = get_cfc_results_all(subjs,data_path); 
%save([data_path 'cfc_results/cfc_all_subj.mat'], 'cfc_all_subj')  
load([data_path 'cfc_results/cfc_all_subj.mat'])
[sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'zstat'); % p value calculations from proportion perms > pixel/num perms
[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);

%https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf
%info for cbrewer color options

%% Fig 3A - Single electrode PLV plot (subj 6 elec 11)

data = cfc_all_subj.s06;
elec_num = 11;
plv_mat = data.norm_plvs{1,find(data.ofc_elecs==elec_num)};

plv_single_elec = figure('Name','zscore PLV one subj one elec')
plv_single_elec.Position = [0, 0, 750, 600]; %best fig size!
contourf(plv_mat);
amp_f_array = 5:5:200;
phase_f_array = 2:20;
% Set locations of ticks and their labels.
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',22,...
'TickDir','in','FontName','Lucida Grande');
title({'\rmPhase-Amplitude Coupling Within','Single OFC Electrode'},'FontSize',30,'FontName','Lucida Grande');
xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Lucida Grande');
ylabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Lucida Grande');

cticks = [-1,1,3,5,7,9,11,13];
clabels = {'-1','1','3','5','7','9','11','13'};
h = colorbar;
h.Ticks = cticks;
h.TickLabels = clabels;
set(get(h, 'ylabel'), 'string', '\fontsize{26}PLV Z-Score','FontSize',22,'FontName','Lucida Grande');
set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path 'single_elec_plv.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

% saveas(plv_single_elec,[fig_path, num2str(elec_num), '_single_elec_plv.pdf'])



%% Fig 3B - Average PLV significant electrodes

sig_plv_group = figure('Name','mean zscore PLV significant elecs across subj')
sig_plv_group.Position = [0, 0, 750, 600]; %best fig size!
contourf(sig_mean_plv);
amp_f_array = 5:5:200;
phase_f_array = 2:20;
% Set locations of ticks and their labels.
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',22,...
'TickDir','in','FontName','Lucida Grande');
title({'\rmOFC Phase-Amplitude Coupling','Across Subjects'},'FontSize',30,'FontName','Lucida Grande');
% title({'\rmStrong Coupling Between Lower Frequency','Phases & Higher Frequency Amplitudes', 'Consistent Across Subjects'},'FontSize',30,'FontName','Lucida Grande');
xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Lucida Grande');
ylabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Lucida Grande');
h = colorbar;
% h.Ticks = cticks;
% h.TickLabels = clabels;
set(get(h, 'ylabel'), 'string', '\fontsize{26}PLV Z-Score','FontSize',22,'FontName','Lucida Grande');
set(gcf,'PaperOrientation','landscape');
% print(gcf,[fig_path 'group_mean_plv_sig_elecs_short_title.pdf'],'-dpdf','-bestfit')

%% Fig 3C - Bar plot PLV values by frequency (block method)
%uses get_plv_freq_info_block_method funct
    %pre-define blocks of interest - 
        %phase_freqs = delta-theta[2 8] alpha[8 13] beta[13 30]
        %amp_freqs = gamma[30 60] hga[60 200] broadgamma[30 200]

%get plot data
        
[all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
    get_plv_freq_block_method(cfc_all_subj, subjs);

freq_specific_plvs = {all_elec_plv_info.delta_theta_gamma_means, all_elec_plv_info.delta_theta_hga_means,...
    all_elec_plv_info.alpha_gamma_means, all_elec_plv_info.alpha_hga_means,...
    all_elec_plv_info.beta_gamma_means,all_elec_plv_info.beta_hga_means};
hga_letters = "HG";
% %define plot colors for phase freq groups - converted to rbg values for plot      
% delta_theta_color = '#de2d26'; % delta_theta_color = '#e34a33'; % 
% alpha_color = '#fc9272'; % alpha_color = '#fdbb84'; % 
% beta_color = '#fee0d2'; % beta_color = '#fee8c8'; % 

%from https://colorbrewer2.org/#type=sequential&scheme=OrRd&n=3 'OrRd' or 'Reds'

%%

plv_by_freq = figure('Name','Mean PLV by Freq');
plv_by_freq.Position = [0,0,750, 600];% plv_by_freq.Position(3:4) = [1000 1000];

x_idxs = [1 2 3];
bar_y_vals = [mean(all_elec_plv_info.delta_theta_gamma_means), mean(all_elec_plv_info.delta_theta_hga_means);...
    mean(all_elec_plv_info.alpha_gamma_means), mean(all_elec_plv_info.alpha_hga_means);...
    mean(all_elec_plv_info.beta_gamma_means),mean(all_elec_plv_info.beta_hga_means)];
b = bar(x_idxs,bar_y_vals,'grouped','FaceColor','flat','EdgeColor', 'w');
set(gca,'XTickLabels', {'\fontname{Lucida Grande}Delta-Theta \bf\delta\theta\rm',...
    '\fontname{Lucida Grande}Alpha \bf\alpha\rm','\fontname{Lucida Grande}Beta \bf\beta\rm'},...
    'FontSize',22,'FontName','Lucida Grande')%...

title({'\rm\fontname{Lucida Grande}Low Frequency Phase Coupling to Gamma \fontsize{34}(\bf\gamma\rm)', '\fontsize{30}& High Gamma (HG) Amplitudes'},...
    'VerticalAlignment','bottom','FontSize',30)
xlabel('\fontname{Lucida Grande}Frequency for Phase (Hz)','FontSize',26)
ylabel('\fontname{Lucida Grande}Phase-Amplitude Coupling (PLV)','FontSize',26)
hold on 
%manually set color data for each bar
b(1,1).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
b(1,1).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
b(1,1).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;[252,174,145]./255;0.9961    0.8980    0.8510

b(1,2).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
b(1,2).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
b(1,2).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;0.9957,0.8784,0.8235

%add amp frequency label to top of bars
xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
ytips1 = b(1,1).YEndPoints; %1x3 Y end points of gamma vals
labels1 = {'\fontsize{20}\fontname{Lucida Grande}\bf\delta\theta\cdot\fontsize{22}\gamma',...
    '\fontsize{20}\fontname{Lucida Grande}\bf\alpha\cdot\fontsize{22}\gamma',...
    '\fontsize{20}\fontname{Lucida Grande}\bf\beta\cdot\fontsize{22}\gamma'};%should be same length as xtips (1x3)
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(1,2).XEndPoints; %1x3 X end points of hga vals
ytips2 = b(1,2).YEndPoints; %1x3 Y end points of hga vals
labels2 = {join(['\fontname{Lucida Grande}\fontsize{20}\bf\delta\theta\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2),...
    join(['\fontname{Lucida Grande}\fontsize{20}\bf\alpha\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2),...
    join(['\fontname{Lucida Grande}\fontsize{20}\bf\beta\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2)};%should be same length as xtips (1x3)
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

% title({'\rm\fontname{Lucida Grande}Low Frequency Phase Coupling with High Frequency Amplitudes:'; 'Gamma (\bf\gamma\rm) Oscillations & High Frequency Activity (\itHFA\rm)'},...
%     'VerticalAlignment','bottom','FontSize',24)
% 
% set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% % set(gcf,'PaperPosition', [0 0 20 10]); %set position of fig to be consistent with pdf sizing
% print(gcf,[fig_path 'barplot.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

% alternative saving options:
% exportgraphics(gcf,[fig_path 'barplot_plv_by_freq.pdf'],"ContentType","vector")
% saveas(plv_by_freq,[fig_path 'barplot_test.pdf'])

% set(gca,'XTickLabels', {'\fontsize{20}\fontname{Lucida Grande}\bf\delta\theta\cdot\fontsize{22}\gamma',...
%     '\fontsize{20}\fontname{Lucida Grande}\bf\alpha\cdot\fontsize{22}\gamma',...
%     '\fontsize{20}\fontname{Lucida Grande}\bf\beta\cdot\fontsize{22}\gamma'},...
%     'FontSize',22,'FontName','Lucida Grande')%...
%% BROADBAND GAMMA ONLY 

plv_by_freq = figure('Name','Mean PLV by Freq');
plv_by_freq.Position = [0,0,700, 600];% plv_by_freq.Position(3:4) = [1000 1000];

x_idxs = [1 2 3];
bar_y_vals = [mean(all_elec_plv_info.delta_theta_broadgamma_means);...
    mean(all_elec_plv_info.alpha_broadgamma_means);...
    mean(all_elec_plv_info.beta_broadgamma_means)];

b = bar(x_idxs,bar_y_vals,'FaceColor','flat','EdgeColor', 'w');

set(gca,'XTickLabels', {'\fontname{Lucida Grande}Delta-Theta',...
    '\fontname{Lucida Grande}Alpha',...
    '\fontname{Lucida Grande}Beta'},...
    'FontSize',22,'FontName','Lucida Grande')

title({'\rm\fontname{Lucida Grande}Coupling Between Low Frequency Phase and'; 'Broadband Gamma Amplitude'},...
    'VerticalAlignment','bottom','FontSize',24)


xlabel('\fontname{Lucida Grande}Frequency for Phase (Hz)','FontSize',24)
ylabel('\fontname{Lucida Grande}Phase-Amplitude Coupling (PLV)','FontSize',24)

hold on 


%manually set color data for each bar
b(1,1).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
b(1,1).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
b(1,1).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;

hold on 
%add sem bars to amp frequency label to top of bars
%calculate sem for bar
b1_data_pts = all_elec_plv_info.delta_theta_broadgamma_means;
b1_sem =  std(all_elec_plv_info.delta_theta_broadgamma_means)/sqrt(15);
b1_upper = mean(b1_data_pts)+b1_sem;
b1_lower = mean(b1_data_pts)-b1_sem;

b2_data_pts = all_elec_plv_info.alpha_broadgamma_means;
b2_sem =  std(all_elec_plv_info.alpha_broadgamma_means)/sqrt(15);
b2_upper = mean(b2_data_pts)+b2_sem;
b2_lower = mean(b2_data_pts)-b2_sem;

b3_data_pts = all_elec_plv_info.beta_broadgamma_means;
b3_sem =  std(all_elec_plv_info.beta_broadgamma_means)/sqrt(15);
b3_upper = mean(b3_data_pts)+b3_sem;
b3_lower = mean(b3_data_pts)-b3_sem;


% err_data = [mean(all_elec_plv_info.delta_theta_broadgamma_means),...
%     mean(all_elec_plv_info.alpha_broadgamma_means),...
%     mean(all_elec_plv_info.beta_broadgamma_means)];
err_data = bar_y_vals;
err_sem = [b1_sem;b2_sem;b3_sem];
% [ngroups, nbars] = size(bar_y_vals);
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);

%get x coords of bars
xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
errorbar(xtips1',err_data,err_sem,'k','linestyle','none','LineWidth', 1.5);%'Color', [.25 .25 .25]);
% err_lower = [b1_lower,b2_lower,b3_lower];
% err_upper = [b1_upper,b2_upper,b3_upper];


% xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
% ytips1 = b(1,1).YEndPoints; %1x3 Y end points of gamma vals


% labels1 = {'\fontname{Lucida Grande}\bf\delta\theta\cdot\gamma',...
%     '\fontname{Lucida Grande}\bf\alpha\cdot\gamma',...
%     '\fontname{Lucida Grande}\bf\beta\cdot\gamma'};%should be same length as xtips (1x3)
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',18)
% set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% % set(gcf,'PaperPosition', [0 0 20 10]); %set position of fig to be consistent with pdf sizing
% print(gcf,[fig_path 'barplot_broadgamma_sem.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page


%% SIG ELEC ONLY BARPLOTS - GAMMA/HGA SPLIT

sig_freq_specific_plvs = {sig_elec_plv_info.delta_theta_gamma_means, sig_elec_plv_info.delta_theta_hga_means,...
    sig_elec_plv_info.alpha_gamma_means, sig_elec_plv_info.alpha_hga_means,...
    sig_elec_plv_info.beta_gamma_means,sig_elec_plv_info.beta_hga_means};
hga_letters = "HG";


%%%%%%split by gamma/high gamma
sig_plv_by_freq = figure('Name','Mean PLV by Freq');
sig_plv_by_freq.Position = [0,0,750, 600];% plv_by_freq.Position(3:4) = [1000 1000];

x_idxs = [1 2 3];
sig_bar_y_vals = [mean(sig_elec_plv_info.delta_theta_gamma_means,'omitnan'), mean(sig_elec_plv_info.delta_theta_hga_means,'omitnan');...
    mean(sig_elec_plv_info.alpha_gamma_means,'omitnan'), mean(sig_elec_plv_info.alpha_hga_means,'omitnan');...
    mean(sig_elec_plv_info.beta_gamma_means,'omitnan'),mean(sig_elec_plv_info.beta_hga_means,'omitnan')];
b = bar(x_idxs,sig_bar_y_vals,'grouped','FaceColor','flat','EdgeColor', 'w');
set(gca,'XTickLabels', {'\fontname{Lucida Grande}Delta-Theta \bf\delta\theta\rm',...
    '\fontname{Lucida Grande}Alpha \bf\alpha\rm','\fontname{Lucida Grande}Beta \bf\beta\rm'},...
    'FontSize',22,'FontName','Lucida Grande')%...

title({'\rm\fontname{Lucida Grande}Sig Elec Low Frequency Phase Coupling to Gamma \fontsize{34}(\bf\gamma\rm)', '\fontsize{30}& High Gamma (HG) Amplitudes'},...
    'VerticalAlignment','bottom','FontSize',30)
xlabel('\fontname{Lucida Grande}Frequency for Phase (Hz)','FontSize',26)
ylabel('\fontname{Lucida Grande}Phase-Amplitude Coupling (PLV)','FontSize',26)
hold on 
%manually set color data for each bar
b(1,1).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
b(1,1).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
b(1,1).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;[252,174,145]./255;0.9961    0.8980    0.8510

b(1,2).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
b(1,2).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
b(1,2).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;0.9957,0.8784,0.8235

%add amp frequency label to top of bars
xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
ytips1 = b(1,1).YEndPoints; %1x3 Y end points of gamma vals
labels1 = {'\fontsize{20}\fontname{Lucida Grande}\bf\delta\theta\cdot\fontsize{22}\gamma',...
    '\fontsize{20}\fontname{Lucida Grande}\bf\alpha\cdot\fontsize{22}\gamma',...
    '\fontsize{20}\fontname{Lucida Grande}\bf\beta\cdot\fontsize{22}\gamma'};%should be same length as xtips (1x3)
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(1,2).XEndPoints; %1x3 X end points of hga vals
ytips2 = b(1,2).YEndPoints; %1x3 Y end points of hga vals
labels2 = {join(['\fontname{Lucida Grande}\fontsize{20}\bf\delta\theta\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2),...
    join(['\fontname{Lucida Grande}\fontsize{20}\bf\alpha\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2),...
    join(['\fontname{Lucida Grande}\fontsize{20}\bf\beta\cdot\rm\fontsize{16}\fontname{Lucida Grande}\it',hga_letters],'',2)};%should be same length as xtips (1x3)
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')


%% SIG ELEC ONLY BARPLOT - BROADGAMMA
%%%%%   BROADGAMMA ONLY

plv_by_freq = figure('Name','Mean PLV by Freq');
plv_by_freq.Position = [0,0,700, 600];% plv_by_freq.Position(3:4) = [1000 1000];

x_idxs = [1 2 3];
bar_y_vals = [mean(sig_elec_plv_info.delta_theta_broadgamma_means,'omitnan');...
    mean(sig_elec_plv_info.alpha_broadgamma_means,'omitnan');...
    mean(sig_elec_plv_info.beta_broadgamma_means,'omitnan')];

b = bar(x_idxs,bar_y_vals,'FaceColor','flat','EdgeColor', 'w');

set(gca,'XTickLabels', {'\fontname{Lucida Grande}Delta-Theta',...
    '\fontname{Lucida Grande}Alpha',...
    '\fontname{Lucida Grande}Beta'},...
    'FontSize',22,'FontName','Lucida Grande')

title({'\rm\fontname{Lucida Grande}Sig Elec Coupling Between Low Frequency Phase and'; 'Broadband Gamma Amplitude'},...
    'VerticalAlignment','bottom','FontSize',24)


xlabel('\fontname{Lucida Grande}Frequency for Phase (Hz)','FontSize',24)
ylabel('\fontname{Lucida Grande}Phase-Amplitude Coupling (PLV)','FontSize',24)

hold on 


%manually set color data for each bar
b(1,1).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
b(1,1).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
b(1,1).CData(3,:) = [0.99,0.8, 0.7];%[254,224,210]./255;

hold on 
%add sem bars to amp frequency label to top of bars
%calculate sem for bar
b1_data_pts = sig_elec_plv_info.delta_theta_broadgamma_means;
b1_sem =  std(sig_elec_plv_info.delta_theta_broadgamma_means,'omitnan')/sqrt(15);
b1_upper = mean(b1_data_pts,'omitnan')+b1_sem;
b1_lower = mean(b1_data_pts,'omitnan')-b1_sem;

b2_data_pts = sig_elec_plv_info.alpha_broadgamma_means;
b2_sem =  std(sig_elec_plv_info.alpha_broadgamma_means,'omitnan')/sqrt(15);
b2_upper = mean(b2_data_pts,'omitnan')+b2_sem;
b2_lower = mean(b2_data_pts,'omitnan')-b2_sem;

b3_data_pts = sig_elec_plv_info.beta_broadgamma_means;
b3_sem =  std(sig_elec_plv_info.beta_broadgamma_means,'omitnan')/sqrt(15);
b3_upper = mean(b3_data_pts,'omitnan')+b3_sem;
b3_lower = mean(b3_data_pts,'omitnan')-b3_sem;


% err_data = [mean(all_elec_plv_info.delta_theta_broadgamma_means),...
%     mean(all_elec_plv_info.alpha_broadgamma_means),...
%     mean(all_elec_plv_info.beta_broadgamma_means)];
err_data = bar_y_vals;
err_sem = [b1_sem;b2_sem;b3_sem];
% [ngroups, nbars] = size(bar_y_vals);
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);

%get x coords of bars
xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
errorbar(xtips1',err_data,err_sem,'k','linestyle','none','LineWidth', 1.5);%'Color', [.25 .25 .25]);




%% Fig 3D - Risk correlation scatter plot DELTA THETA GAMMA RISK TRENDING BUT NOT ANYTHIGN FOR NOW
% mean_deltatheta_gamma = all_elec_plv_info.delta_theta_hga_means;
mean_deltatheta_gamma = sig_elec_plv_info.delta_theta_broadgamma_means(~isnan(sig_elec_plv_info.delta_theta_broadgamma_means)); %SIG
mean_deltatheta_gamma = log(mean_deltatheta_gamma);
% nonsig_subj_idx = find(isnan(mean_deltatheta_gamma));
risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
risk_prefs = risk_prefs(~isnan(mean_deltatheta_gamma)); %remove nonsig subj 


lm = fitlm(mean_deltatheta_gamma,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(mean_deltatheta_gamma.',risk_prefs.','type','pearson');

scatter(mean_deltatheta_gamma,risk_prefs)
xlabel('log mean plv btwn delta-theta and hga')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))

%delta theta gamma p 0.366, hga p 0.22 broadgamma p 0.249
%alpha p 0.554 beta p 0.419
% 
% %% Fig 3E - Coherence correlation scatter plot
% coh_header = {'RiskPref', 'coh_delta_theta', 'coh_delta',	'coh_theta',	'coh_alpha',  'coh_beta'};
% coh_delta_theta = mean_coh(:,3)'; %load by hand for now...
% mean_delta_theta_gamma_plvs = all_elec_plv_info.delta_theta_hga_means; %0.363
% 
% lm = fitlm(mean_delta_theta_gamma_plvs,coh_delta_theta)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(mean_delta_theta_gamma_plvs.',coh_delta_theta.','type','pearson');
% 
% 
% scatter(mean_delta_theta_gamma_plvs,coh_delta_theta)
% ylabel('Mean Delta-Theta Coherence')
% xlabel('Mean Delta-Theta Gamma PLV')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))

% %% CHECK FOR CORRS WITH SIG SUBJ ONLY - nothing for risk at all
% sig_subj = subjs(1,find([sig_elecs_by_subj{1,:}]~=0));
% sig_subj_idx = find([sig_elecs_by_subj{1,:}]~=0);
% 
% sig_subj_risks = risk_prefs(sig_subj_idx);
% 
% coh_delta_theta = mean_coh(:,4)'; %load by hand for now... 0.359 gamma delta theta
% sig_subj_cohs = coh_delta_theta(sig_subj_idx);
% 
% mean_deltatheta_gamma = all_elec_plv_info.delta_theta_broadgamma_means;
% sig_subj_plvs = mean_deltatheta_gamma(sig_subj_idx);
% 
% lm = fitlm(sig_subj_plvs,sig_subj_cohs)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(sig_subj_plvs.',sig_subj_cohs.','type','pearson');
% 
% scatter(sig_subj_plvs,sig_subj_cohs)
% xlabel('sig_subj_plvs')
% ylabel('sig_subj_cohs')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))


%% fds
% %% Fig 3C - bar plot plv values by freq (cluster method)
% %uses get_plv_freq_info_clust_method funct
% 
% %%%% extract the frequencies of amp/phase that are in the most significant
% %%%% contours
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
% 
% 
% [all_elec_plv_info,sig_elec_plv_info] = get_plv_freq_info_clust_method(cfc_all_subj, subjs, sig_amp_idx,sig_phase_idx);
% 
% 
% %%%%%% barplot of all elecs
% bar_x = [1, 2, 3, 4, 5];
% %bar1 = total plv 
% bar1_x = bar_x(1);
% bar1_y = mean(all_elec_plv_info.mean_plv_norm);
% scat1_y = all_elec_plv_info.mean_plv_norm;
% scat1_x = ones(1,length(scat1_y));
% %bar2 = all sig clust 
% bar2_x = bar_x(2);
% bar2_y = mean(all_elec_plv_info.mean_total_cluster);
% scat2_y = all_elec_plv_info.mean_total_cluster;
% scat2_x = ones(1,length(scat2_y))*bar2_x;
% % %bar3 = beta
% % bar3_x = bar_x(3);
% % bar3_y = mean(all_elec_plv_info.mean_deltatheta_beta);
% % scat3_y = all_elec_plv_info.mean_deltatheta_beta;
% % scat3_x = ones(1,length(scat3_y))*bar3_x;
% %bar3 = gamma
% bar3_x = bar_x(3);
% bar3_y = mean(all_elec_plv_info.mean_deltatheta_gamma);
% scat3_y = all_elec_plv_info.mean_deltatheta_gamma;
% scat3_x = ones(1,length(scat3_y))*bar3_x;
% %bar4 = hga
% bar4_x = bar_x(4);
% bar4_y = mean(all_elec_plv_info.mean_deltatheta_hga);
% scat4_y = all_elec_plv_info.mean_deltatheta_hga;
% scat4_x = ones(1,length(scat4_y))*bar4_x;
% %bar5 = broadgamma
% bar5_x = bar_x(5);
% bar5_y = mean(all_elec_plv_info.mean_deltatheta_broadgamma);
% scat5_y = all_elec_plv_info.mean_deltatheta_broadgamma;
% scat5_x = ones(1,length(scat5_y))*bar5_x;
% 
% %%%%%%bar and points
% plv_by_freq = figure('Name','All Subj Mean PLV by Freq');
% % ylim([0,2.5])
% bar(bar1_x,bar1_y)
% hold on
% scatter(scat1_x,scat1_y)
% hold on 
% bar(bar2_x,bar2_y)
% hold on 
% scatter(scat2_x,scat2_y)
% hold on
% bar(bar3_x,bar3_y)
% hold on 
% scatter(scat3_x,scat3_y)
% hold on
% bar(bar4_x,bar4_y)
% hold on 
% scatter(scat4_x,scat4_y)
% hold on
% bar(bar5_x,bar5_y)
% hold on
% scatter(scat5_x,scat5_y)
% 
% %% bar only plot
% %%%%%%%%%%baronly 
% x_labels = {'All Freqs (20-200Hz)','Cluster Freqs (5-170Hz)','Gamma(30-70Hz)','HGA (70-200Hz)','Broadband Gamma(30-200Hz)'};
% plv_by_freq = figure('Name','bar only All Subj Mean PLV by Freq');
% title('Delta-Theta Phase Locking Across Frequencies')
% % ylim([0,2.5])
% bar(bar1_x,bar1_y,'FaceColor',"#0072BD")
% % hold on
% % scatter(scat1_x,scat1_y)
% hold on 
% bar(bar2_x,bar2_y,'FaceColor',green)%"#7E2F8E"
% % hold on 
% % scatter(scat2_x,scat2_y)
% hold on
% bar(bar3_x,bar3_y,'FaceColor',"#EDB120") %"#D95319"
% % hold on 
% % scatter(scat3_x,scat3_y)
% hold on
% bar(bar4_x,bar4_y,'FaceColor',red) %"#A2142F"
% % hold on 
% % scatter(scat4_x,scat4_y)
% hold on
% bar(bar5_x,bar5_y,'FaceColor',brown)
% ax=gca;
% ax.XTick = bar_x;
% ax.XTickLabels=x_labels;
% ax.XTickLabelRotation = 45;
% % hold on
% % scatter(scat5_x,scat5_y)
% 
% %stderror = std( data ) / sqrt( length ) for baronly
% 
% [222,45,38]./255;
% % red = [204 37 41]./255;%[.204 .37 .41]
% % % red = [211 94 96]./255;
% % green = [62 150 81]./255;
% % brown = [146 36 40]./255;
% 

%% broadband only 
% 
% plv_by_freq = figure('Name','Mean PLV by Freq');
% x_idxs = [1 2 3];
% bar_y_vals = [mean(all_elec_plv_info.delta_theta_broadgamma_means);...
%     mean(all_elec_plv_info.alpha_broadgamma_means);...
%     mean(all_elec_plv_info.beta_broadgamma_means)];
% 
% b = bar(x_idxs,bar_y_vals,'FaceColor','flat','EdgeColor', 'w');
% set(gca, 'XTickLabels', {'Delta-Theta','Alpha','Beta'})
% xlabel('Frequency of Phase (Hz)')
% ylabel('Mean PLV')
% hold on 
% %manually set color data for each bar
% b(1,1).CData(1,:) = [0.8706,0.1765,0.1490];%[222,45,38]./255;
% b(1,1).CData(2,:) = [0.9882,0.5725,0.4471];%[252,146,114]./255;
% b(1,1).CData(3,:) = [0.9961,0.8784,0.8235];%[254,224,210]./255;
% 
% 
% %add amp frequency label to top of bars
% xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
% ytips1 = b(1,1).YEndPoints; %1x3 Y end points of gamma vals
% labels1 = {'\bf\delta\theta\rm-\bf\gamma','\bf\alpha\rm-\bf\gamma','\bf\beta\rm-\bf\gamma'};%should be same length as xtips (1x3)
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',18)


%% Supplementary Figures 

% Histogram of permutation values 


surr = squeeze(cfc_all_subj.s02.perm_plvs{1,1}(1,1,:));

plv_zscore = cfc_all_subj.s02.norm_real_plvs{1,1}(1,1);

perm_fig = figure('Name','Single Electrode Single Pixel Permutation Distribution')
histogram(surr,20,'FaceColor','#D95319')
hold on 
xline(plv_zscore,'--k','Observed PLV zscore','LineWidth',1.5)
ylabel('Frequency')
xlabel('Surrogate Phase Locking Values')
title('Example Significant Permutation Distribution (n=1000)')
saveas(perm_fig,[fig_path 'permutation_example.pdf'])

% Another way to calculate pvalue - gives slightly different result...
% F = @(x)(exp (-0.5*(x.^2))./sqrt (2*pi));
% p = quad(F, plv_zscore, 100);

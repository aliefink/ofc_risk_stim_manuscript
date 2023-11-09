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
load([data_path 'cfc_results/cfc_all_subj_sig_info.mat'])
load([data_path 'cfc_results/sig_elecs_by_subj.mat'])
load([data_path 'cfc_results/nonsig_elecs_by_subj.mat'])
load([data_path 'cfc_results/find_sig_elecs_header.mat'])

% [sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'prop'); % p value calculations from proportion perms > pixel/num perms
[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);

%https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf
%info for cbrewer color options


%% Extract mean plv data by freq band


[all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
    get_plv_freq_block_method(cfc_all_subj, subjs);

%% Fig 2A - Single electrode PLV plot (subj 8 elec 94)
subj_id = 's06';
elec_id = '11';
data = cfc_all_subj.s06; % s06 11; s08 94; s04 80
elec_num = 11;
plv_mat = data.norm_plvs{1,find(data.ofc_elecs==elec_num)};

plv_single_elec = figure('Name','zscore PLV one subj one elec')
plv_single_elec.Position = [0, 0, 750, 600]; %best fig size!
contourf(plv_mat,5);
amp_f_array = 5:5:200;
phase_f_array = 2:20;
% Set locations of ticks and their labels.
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',22,...
'TickDir','in','FontName','Lucida Grande');
title({'\rmPhase-Amplitude Coupling Within','Single OFC Electrode'},'FontSize',30,'FontName','Lucida Grande');
xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Lucida Grande');
ylabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Lucida Grande');

% cticks = [-1,1,3,5,7,9,11,13];
% clabels = {'-1','1','3','5','7','9','11','13'};
h = colorbar;
% h.Ticks = cticks;
% h.TickLabels = clabels;
set(get(h, 'ylabel'), 'string', '\fontsize{26}PLV Z-Score','FontSize',22,'FontName','Lucida Grande');
set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path subj_id '_' elec_id '_single_elec_plv.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

% saveas(plv_single_elec,[fig_path, num2str(elec_num), '_single_elec_plv.pdf'])



%% Fig 2B - Average PLV significant electrodes

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

%% Fig 2C - Line plot of frequency distributions

%phase frequency distribution by mean plvs
plv_distr_phase_freqs = mean(sig_mean_plv,1); %mean along columns - phase means

figure
plot(phase_f_array,plv_distr_phase_freqs)

%amp frequency distribution by mean plvs
figure
plv_distr_amp_freqs = mean(sig_mean_plv,2)'; %mean along rows - amp means
plot(amp_f_array,plv_distr_amp_freqs) %max peak @ 35 second @ 90


%to do with just broadgamma - 
broadgamma = [30 200];
broadgamma_idx = find((amp_f_array>=broadgamma(1)&amp_f_array<=broadgamma(2)));
broadgamma_plvs_mean = sig_mean_plv(broadgamma_idx,:);

broadgamma_distr_phase_freqs = mean(broadgamma_plvs_mean,1);

plot(phase_f_array,broadgamma_distr_phase_freqs)


%% Fig 2D - Bar plot PLV values by frequency (block method)
%uses get_plv_freq_info_block_method funct
    %pre-define blocks of interest - 
        %phase_freqs = delta-theta[2 8] alpha[8 13] beta[13 30]
        %amp_freqs = gamma[30 60] hga[60 200] broadgamma[30 200]


% %define plot colors for phase freq groups - converted to rbg values for plot      
% delta_theta_color = '#de2d26'; % delta_theta_color = '#e34a33'; % 
% alpha_color = '#fc9272'; % alpha_color = '#fdbb84'; % 
% beta_color = '#fee0d2'; % beta_color = '#fee8c8'; % 

%from https://colorbrewer2.org/#type=sequential&scheme=OrRd&n=3 'OrRd' or 'Reds'


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



%% SUPPLEMENTARY FIGURES

%% Supplementary Figures 

% Histogram of permutation values 

subj_id = 's06';
elec_id = '11';
data = cfc_all_subj.s06; % s06 11; s08 94; s04 80
elec_num = 11;
elec_idx = find(data.ofc_elecs==elec_num);
surr = squeeze(cfc_all_subj.s06.perm_plvs{1,elec_idx}(12,11,:));

plv_pixel = cfc_all_subj.s06.real_plvs{1,elec_idx}(12,11);

perm_fig = figure('Name','Single Electrode Single Pixel Permutation Distribution')
histogram(surr,20,'FaceColor','#D95319')
hold on 
xline(plv_pixel,'--k','Observed PLV','LineWidth',1.5)
ylabel('Count')
xlabel('Surrogate Phase Locking Values')
title('Example Significant Permutation Distribution (n=1000)')
% saveas(perm_fig,[fig_path 'SOM_methods/sig_permutation_example.pdf'])

%%%%NOTSIG DISTR
surr = squeeze(cfc_all_subj.s06.perm_plvs{1,elec_idx}(12,15,:));

plv_pixel = cfc_all_subj.s06.real_plvs{1,elec_idx}(12,15);

perm_fig = figure('Name','Single Electrode Single Pixel Permutation Distribution')
histogram(surr,20,'FaceColor','#D95319')
hold on 
xline(plv_pixel,'--k','Observed PLV','LineWidth',1.5)
ylabel('Count')
xlabel('Surrogate Phase Locking Values')
title('Example Not Significant Permutation Distribution (n=1000)')
% saveas(perm_fig,[fig_path 'SOM_methods/nonsig_permutation_example.pdf'])
%% Risk correlation scatter plot DELTA THETA GAMMA RISK TRENDING BUT NOT ANYTHIGN FOR NOW
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



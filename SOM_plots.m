
%% SUPPLEMENTARY FIGURES

% CFC SOM Fig 1:
    %nonsig PLV comodulogram
    %broadgamma bar plot with phase freqs
    %PLV risk correlation scatter plot
    

% CFC SOM Fig 2:
    %sig surrogate histogram 
    %sig fdr correction
    %nonsig surrogate histogram
    %nonsig fdr correction
    %cluster graphic if needed later...

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
    
%subj_ids in cfc dataset 
subjs = {'s02','s04','s06','s07','s08','s09','s10','s12','s13','s14','s15','s16','s17','s18','s19'};
load([data_path 'cfc_results/cfc_all_subj.mat']) 
[sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'prop','no'); % can select pval 'zstat' or 'prop'

%% CFC SOM Fig 1
%% SOM Fig 1a NONSIG CFC COMODULOGRAM

[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);

%%%%NONSIG

nonsig_plv_group = figure('Name','Not Significant Elec PLVs')
nonsig_plv_group.Position = [0, 0, 750, 600]; %best fig size!

contourf(nonsig_mean_plv);

amp_f_array = 5:5:200;
phase_f_array = 2:20;


% Set locations of ticks and their labels.
lim = max(sig_mean_plv,[],'all');

set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',22,...
'TickLength',[0 0],'FontName','Arial', 'CLim', [0 lim] );

xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Arial');
ylabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Arial');


h = colorbar;
set(get(h, 'ylabel'), 'string', '\fontsize{26}PLV Z-Score','FontSize',22,'FontName','Arial');
set(gca, 'box', 'off')
set(gcf,'PaperOrientation','landscape'); 


% print(gcf,[fig_path 'nonsig_group_plv_comod.pdf'],'-dpdf','-bestfit')

% %%%%%summary stats for fig:
% mean(nonsig_mean_plv,'all')
% std(nonsig_mean_plv,[],'all','omitnan')
% max(nonsig_mean_plv,[],'all')
% [row,col] = find(nonsig_mean_plv == max(nonsig_mean_plv,[],'all')) %row = amp idx, col = phase idx
% phase_f_array(col)
% amp_f_array(row)

%% SOM Fig 1b BAR PLOT BROADGAMMA X PHASE FREQ BANDS
    %Bar plot PLV values by frequency (block method)
%uses get_plv_freq_info_block_method funct
    %pre-define blocks of interest - 
        %phase_freqs = delta-theta[2 8] alpha[8 13] beta[13 30]
        %amp_freqs = gamma[30 70] hga[70 200] broadgamma[30 200]



plv_by_freq = figure('Name','Mean PLV by Freq');
plv_by_freq.Position = [0,0,700, 600];% plv_by_freq.Position(3:4) = [1000 1000];

x_idxs = [1 2 3];
bar_y_vals = [mean(sig_elec_plv_info.delta_theta_broadgamma_means,'omitnan');...
    mean(sig_elec_plv_info.alpha_broadgamma_means,'omitnan');...
    mean(sig_elec_plv_info.beta_broadgamma_means,'omitnan')];

b = bar(x_idxs,bar_y_vals,'FaceColor','flat','EdgeColor', 'w');

set(gca,'FontSize',18,'FontName','Arial','TickLength',[0 0])
xticklabels({'\fontsize{28}\fontname{Arial}\bf\delta\theta\cdot\gamma',...
    '\fontsize{28}\fontname{Arial}\bf\alpha\cdot\gamma',...
    '\fontsize{28}\fontname{Arial}\bf\beta\cdot\gamma'})
% xticklabels({'\fontsize{24}Delta-Theta',...
%     '\fontsize{24}Alpha',...
%     '\fontsize{24}Beta'})

xlabel('\fontname{Arial}Phase-Amplitude Frequencies (Hz)','FontSize',26)
ylabel('\fontname{Arial}PLV Z-Score','FontSize',26)

hold on 


%manually set color data for each bar
%saez lab colors: 255 172 77; 254 189 113; 255 204 147;  255,222,184; 255,238,219
b(1,1).CData(1,:) = [255 172 77]./255;%[203,24,29]./255;brightR[203,24,29]deepR[165,15,21]Or[166,54,3]
b(1,1).CData(2,:) = [255 204 147]./255;%[252,146,114]./255;Or[241,105,19]
b(1,1).CData(3,:) = [255 238 219]./255;%[254,224,210]./255;Or[253,174,107]

hold on 
%add sem bars to amp frequency label to top of bars
%calculate sem for bar
b1_data_pts = sig_elec_plv_info.delta_theta_broadgamma_means;
b1_sem =  std(sig_elec_plv_info.delta_theta_broadgamma_means,'omitnan')/sqrt(14);
b1_upper = mean(b1_data_pts,'omitnan')+b1_sem;
b1_lower = mean(b1_data_pts,'omitnan')-b1_sem;

b2_data_pts = sig_elec_plv_info.alpha_broadgamma_means;
b2_sem = std(sig_elec_plv_info.alpha_broadgamma_means,'omitnan')/sqrt(14);
b2_upper = mean(b2_data_pts,'omitnan')+b2_sem;
b2_lower = mean(b2_data_pts,'omitnan')-b2_sem;

b3_data_pts = sig_elec_plv_info.beta_broadgamma_means;
b3_sem =  std(sig_elec_plv_info.beta_broadgamma_means,'omitnan')/sqrt(14);
b3_upper = mean(b3_data_pts,'omitnan')+b3_sem;
b3_lower = mean(b3_data_pts,'omitnan')-b3_sem;

err_data = bar_y_vals;
err_sem = [b1_sem;b2_sem;b3_sem];

%get x coords of bars
xtips1 = b(1,1).XEndPoints; %1x3 X end points of gamma vals
%add errorbars using errorbar function
errorbar(xtips1',err_data,err_sem,'k','linestyle','none','LineWidth', 1.5,'CapSize',5);%'Color', [.25 .25 .25]);
set(gca, 'box', 'off')
set(gcf,'PaperOrientation','landscape');
% print(gcf,[fig_path 'barplot_broadgamma_sem.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page


%%%%significance testing for barplot 
anova_data = [sig_elec_plv_info.delta_theta_broadgamma_means; ...
    sig_elec_plv_info.alpha_broadgamma_means; ...
    sig_elec_plv_info.beta_broadgamma_means]';
[p,tbl,stats] = anova1(anova_data)



%% SOM Fig 1c PLV x RISK CORRELATION
    %Risk correlation scatter plot DELTA THETA GAMMA RISK

mean_deltatheta_gamma = sig_elec_plv_info.delta_theta_broadgamma_means(~isnan(sig_elec_plv_info.delta_theta_broadgamma_means)); %SIG
%mean_deltatheta_gamma = log(mean_deltatheta_gamma); %option to log transform
risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
risk_prefs = risk_prefs(~isnan(mean_deltatheta_gamma)); %remove nonsig subj 
% risk_prefs = log(risk_prefs);


lm = fitlm(mean_deltatheta_gamma,risk_prefs)
rsq = round(lm.Rsquared.Ordinary,2);
p_val = round(table2array(lm.Coefficients(2,4)),3);
r=round(corr(mean_deltatheta_gamma.',risk_prefs.','type','pearson'),2);

ip_deltatheta_gamma_fig = figure('Name','Indifference Point x PLV');
ip_deltatheta_gamma_fig.Position = [0,0,750, 600];% plv_by_freq.Position(3:4) = [1000 1000];

scatter(mean_deltatheta_gamma,risk_prefs,150,[255 172 77]./255,'filled')
hold on
h = lsline;
set(h,'LineWidth',2.5,'Color','k','LineStyle','--'); %[193 155 211]./255
annotation('textbox',[.75 .75 .15 .125],'string',strcat('R^{2} = '," ",string(rsq),...
    '\newline\itr\rm ='," ",string(r), '\newline\itp\rm = '," ", string(p_val)),...
    'FontSize',20,'FontName','Arial','LineWidth',1,'LineStyle','none',...
    'VerticalAlignment','middle')


set(gca,'FontSize',18,'FontName','Arial','TickLength',[0 0])
xlabel('ln PLV Z-Score (\bf\delta\theta\cdot\gamma\rm)','FontSize',24)
% xlabel('PLV Z-Score (Delta-Theta:Gamma)','FontSize',24)

ylabel('Indifference Point','FontSize',24)
set(gcf,'PaperOrientation','landscape');
set(gca, 'box', 'off')

%run below to save:
% print(gcf,[fig_path 'plv_risk_broadgamma.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page


% outlier_idx = find(mean_deltatheta_gamma == max(mean_deltatheta_gamma));
% mean_deltatheta_gamma(outlier_idx) = [];
% risk_prefs(outlier_idx) = [];


%% CFC SOM Fig 2 
    
%fig colors: %1 blue 2 purple 3 yellow 4 lavender 5 dusty rose
saez_cmap = {'#00467e', '#561C46', '#FFAC4D', '#C19BD3','#CDB6B6'}; 

%% SOM Fig 2a SIG PIXEL SURROGATES HISTOGRAM

subj_id = 's06';
elec_id = '11';
data = cfc_all_subj.s06; % s06 11; s08 94; s04 80
elec_num = 11;
elec_idx = find(data.ofc_elecs==elec_num);
surr = squeeze(cfc_all_subj.s06.perm_plvs{1,elec_idx}(12,11,:)); %60hz amp 12 hz phase

plv_pixel = cfc_all_subj.s06.real_plvs{1,elec_idx}(12,11);

perm_fig = figure('Name','Sig Pixel Permutation Distribution')
histogram(surr,20,'FaceColor','#00467e','FaceAlpha',0.75)
set(gca,'FontSize',18,'FontName','Arial','TickLength',[0 0])
hold on 
xline(plv_pixel,'--k',{'Observed', 'PLV'}, 'LineWidth',1.5,'FontSize',20,...
    'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
ylabel('Count','FontSize',24)
xlabel('Surrogate PLV','FontSize',24)
set(gcf,'PaperOrientation','landscape');
set(gca, 'box', 'off')
% saveas(perm_fig,[fig_path 'sig_permutation_example.pdf'])

%%%%%summary stats for fig 
% plv_pixel %real plv 
% (plv_pixel-mean(surr))/std(surr) %zscore plv
% mean(surr) %surrogate mean
% std(surr) %surrogate std
% max(surr) %surrogate max
% sum(surr>plv_pixel) %num surr > real 
% sum(surr>plv_pixel)/760 % pvalue 


%% SOM Fig 2b NONSIG PIXEL SURROGATES HISTOGRAM
%%%%NOTSIG DISTR

surr = squeeze(cfc_all_subj.s06.perm_plvs{1,elec_idx}(12,15,:));

plv_pixel = cfc_all_subj.s06.real_plvs{1,elec_idx}(12,15); %60 hz amp and 16 hz phase

perm_fig = figure('Name','Sig Pixel Permutation Distribution')
histogram(surr,20,'FaceColor','#00467e','FaceAlpha',0.75)
set(gca,'FontSize',18,'FontName','Arial','TickLength',[0 0])
hold on 
xline(plv_pixel,'--k',{'Observed', 'PLV'}, 'LineWidth',1.5,'FontSize',20,...
    'LabelOrientation','horizontal','LabelHorizontalAlignment','right')
ylabel('Count','FontSize',24)
xlabel('Surrogate PLV','FontSize',24)
set(gcf,'PaperOrientation','landscape');
set(gca, 'box', 'off')
% saveas(perm_fig,[fig_path 'nonsig_permutation_example.pdf'])

% %%%%%summary stats for fig 
% plv_pixel %real plv 
% (plv_pixel-mean(surr))/std(surr) %zscore plv
% mean(surr) %surrogate mean
% std(surr) %surrogate std
% max(surr) %surrogate max
% sum(surr>plv_pixel) %num surr > real 
% sum(surr>plv_pixel)/760 % pvalue 


%% SOM FIG 2C - SIG ELECTRODE FDR Correction Schematic 
%fig colors: %1 blue 2 purple 3 yellow 4 lavender 5 dusty rose
saez_cmap = {'#00467e', '#561C46', '#FFAC4D', '#C19BD3','#CDB6B6'}; 

%sig = s06 e 11
subj_id = 's06';
elec_id = '11';
data = cfc_all_subj.s06; % s06 11; s08 94; s04 80
elec_num = 11;

%uncorrected pvalues
og_pvals = data.prop_pvals{1,find(data.ofc_elecs==elec_num)};
og_pval_vec = reshape(og_pvals,1,[]);
num_og_sig = string(sum(og_pval_vec<0.05));
%corrected pvalues
fdr_pvals = sig_elecs_by_subj{6,find(strcmp(subjs,subj_id))}{find(sig_elecs_by_subj{2,find(strcmp(subjs,subj_id))}==elec_num),1}; %adj pval is row 6,subj is col, elec is row in cell 
fdr_pval_vec = reshape(fdr_pvals,1,[]);
num_fdr_sig = string(sum(fdr_pval_vec<0.05));
%total number of pixels
num_pixels = string(length(fdr_pval_vec));

pval_fig = figure('Name','pvalue distributions') 
pval_fig.Position = [0, 0, 750, 600]; %best fig size!
histogram(og_pval_vec,'FaceAlpha',0.8,'FaceColor','#561C46','Normalization','count','BinWidth',0.05);%,'DisplayName','Uncorrected p values');
hold on;
histogram(fdr_pval_vec,'FaceAlpha',0.8,'FaceColor','#FFAC4D','Normalization','count','BinWidth',0.05);%,'DisplayName','FDR Adjusted p values\newline num sig'); %set bin width to alpha for count of < alpha
hold on;
xline(0.05,'--k',...
    '\fontname{Arial}\fontsize{20}Significance Threshold\newline\fontsize{24}\bf\alpha\rm\fontname{Arial}\fontsize{22} = 0.05',...
    'LineWidth',1.5,'Alpha',1,...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',22,'FontName','Arial');
xlabel('PLV \itp\rm value','FontSize',28,'FontName','Arial');
ylabel('Count','FontSize',28,'FontName','Arial');
ylim([0 500])
% xlim([0,max(fdr_pval_vec)])
legend({'Uncorrected \itp\rm values',...
    'FDR Corrected \itp\rm values'},...
    'FontSize',24,'FontName','Arial','Location','northeast')

set(gca, 'box', 'off')

set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  'significant_elec_pval_correction.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

%%%%summary stats for figure
% mean(og_pval_vec,'omitnan')
% std(og_pval_vec,'omitnan')
% num_og_sig
% mean(fdr_pval_vec,'omitnan')
% std(fdr_pval_vec,'omitnan')
% num_fdr_sig

%% SOM FIG2D NONSIG ELECTRODE FDR Correction
%nonsig = s02 e 22
subj_id = 's02';
elec_id = '21';
data = cfc_all_subj.s02; % s06 11; s08 94; s04 80
elec_num = 21;

%uncorrected pvalues
og_pvals = data.prop_pvals{1,find(data.ofc_elecs==elec_num)};
og_pval_vec = reshape(og_pvals,1,[]);
num_og_sig = string(sum(og_pval_vec<0.05));
%corrected pvalues
fdr_pvals = nonsig_elecs_by_subj{6,find(strcmp(subjs,subj_id))}{find(nonsig_elecs_by_subj{2,find(strcmp(subjs,subj_id))}==elec_num),1}; %adj pval is row 6,subj is col, elec is row in cell 
fdr_pval_vec = reshape(fdr_pvals,1,[]);
num_fdr_sig = string(sum(fdr_pval_vec<0.05));
%total number of pixels
num_pixels = string(length(fdr_pval_vec));

pval_fig = figure('Name','pvalue distributions') 
pval_fig.Position = [0, 0, 750, 600]; %best fig size!
histogram(og_pval_vec,'FaceAlpha',0.8,'FaceColor','#561C46','Normalization','count','BinWidth',0.05);%,'DisplayName','Uncorrected p values');
hold on;
histogram(fdr_pval_vec,'FaceAlpha',0.8,'FaceColor','#FFAC4D','Normalization','count','BinWidth',0.05);%,'DisplayName','FDR Adjusted p values\newline num sig'); %set bin width to alpha for count of < alpha
hold on;
xline(0.05,'--k',...
    '\fontname{Arial}\fontsize{20}Significance Threshold\newline\fontsize{24}\bf\alpha\rm\fontname{Arial}\fontsize{22} = 0.05',...
    'LineWidth',1.5,'Alpha',1,...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
xlabel('PLV \itp\rm value','FontSize',28,'FontName','Arial');
ylabel('Count','FontSize',28,'FontName','Arial');
ylim([0 120])
legend({'Uncorrected \itp\rm values',...
    'FDR Corrected \itp\rm values'},...
    'FontSize',24,'FontName','Arial','Location','northeast')
set(gca,'TickLength',[0 0],'FontSize',22,'FontName','Arial');

set(gca, 'box', 'off')

set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  'nonsig_elec_pval_correction.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page


%%%%summary stats for figure
% mean(og_pval_vec,'omitnan')
% std(og_pval_vec,'omitnan')
% num_og_sig
% mean(fdr_pval_vec,'omitnan')
% std(fdr_pval_vec,'omitnan')
% num_fdr_sig


%% MISC ARCHIVE

% pval_fig = figure('Name','pvalue distributions') 
% pval_fig.Position = [0, 0, 750, 600]; %best fig size!
% histogram(og_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05,...
%     'DisplayName',join(['\fontname{Arial}\fontsize{18}Uncorrected \itp\rm values \newline '...
%     num_og_sig '/' num_pixels ' < \fontsize{22}\bf\alpha\rm\fontsize{18}'],''));
% hold on;
% xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
%     '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
%     'LabelOrientation','horizontal')
% set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
% xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
% ylabel('Count','FontSize',22,'FontName','Arial');
% legend()
% set(gca, 'box', 'off')
% set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% % print(gcf,[fig_path  '/SOM_methods/NONsignificant_pval_uncorrected_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page
% 
% 
% 
% pval_fig = figure('Name','pvalue distributions') 
% pval_fig.Position = [0, 0, 750, 600]; %best fig size!
% histogram(fdr_pval_vec,'FaceColor','#D95319','FaceAlpha',0.8,...
%     'Normalization','count','BinWidth',0.05,'DisplayName',...
%     join(['FDR Corrected \itp\rm values \newline' num_fdr_sig ...
%     '/' num_pixels ' < \fontsize{22}\bf\alpha'],'')); %set bin width to alpha for count of < alpha
% hold on;
% xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
%     '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
%     'LabelOrientation','horizontal')
% set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
% xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
% ylabel('Count','FontSize',22,'FontName','Arial');
% % xlim([0,max(fdr_pval_vec)])
% legend()
% 
% set(gca, 'box', 'off')
% 
% set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% % print(gcf,[fig_path  '/SOM_methods/NONsignificant_pval_fdr_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page


% %% Separate the histograms
% %uncorrected pvalues
% og_pvals = data.prop_pvals{1,find(data.ofc_elecs==elec_num)};
% og_pval_vec = reshape(og_pvals,1,[]);
% num_og_sig = string(sum(og_pval_vec<0.05));
% 
% pval_fig = figure('Name','pvalue distributions') 
% pval_fig.Position = [0, 0, 750, 600]; %best fig size!
% histogram(og_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05,...
%     'DisplayName',join(['\fontname{Arial}\fontsize{18}Uncorrected \itp\rm values \newline '...
%     num_og_sig '/' num_pixels ' < \fontsize{22}\bf\alpha\rm\fontsize{18}'],''));
% hold on;
% xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
%     '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
%     'LabelOrientation','horizontal')
% set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
% xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
% ylabel('Count','FontSize',22,'FontName','Arial');
% legend()
% set(gca, 'box', 'off')
% set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% % print(gcf,[fig_path  '/SOM_methods/significant_pval_uncorrected_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page
% 
% 
% 
% %corrected pvalues
% fdr_pvals = sig_elecs_by_subj{6,find(strcmp(subjs,subj_id))}{find(sig_elecs_by_subj{2,find(strcmp(subjs,subj_id))}==elec_num),1}; %adj pval is row 6,subj is col, elec is row in cell 
% fdr_pval_vec = reshape(fdr_pvals,1,[]);
% num_fdr_sig = string(sum(fdr_pval_vec<0.05));
% %total number of pixels
% num_pixels = string(length(fdr_pval_vec));
% 
% pval_fig = figure('Name','pvalue distributions') 
% pval_fig.Position = [0, 0, 750, 600]; %best fig size!
% histogram(fdr_pval_vec,'FaceColor','#D95319','FaceAlpha',0.8,...
%     'Normalization','count','BinWidth',0.05,'DisplayName',...
%     join(['FDR Corrected \itp\rm values \newline' num_fdr_sig ...
%     '/' num_pixels ' < \fontsize{22}\bf\alpha'],'')); %set bin width to alpha for count of < alpha
% hold on;
% xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
%     '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
%     'LabelOrientation','horizontal')
% set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
% xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
% ylabel('Count','FontSize',22,'FontName','Arial');
% % xlim([0,max(fdr_pval_vec)])
% legend()
% 
% set(gca, 'box', 'off')
% 
% set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% % print(gcf,[fig_path  '/SOM_methods/significant_pval_fdr_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page
% 
% 

% 
% 
% %% Fig 3E - Coherence correlation scatter plot
% coh_header = {'RiskPref', 'coh_delta_theta', 'coh_delta',	'coh_theta',	'coh_alpha',  'coh_beta'};
% load('/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/data/coh_data/OrG-OrG-Decision_coh_all_lowfreq.mat')
% coh_delta_theta = mean_coh(:,3)';
% % mean_delta_theta_gamma_plvs = all_elec_plv_info.delta_theta_hga_means; %0.363
% mean_deltatheta_gamma = sig_elec_plv_info.delta_theta_broadgamma_means(~isnan(sig_elec_plv_info.delta_theta_broadgamma_means)); %SIG
% mean_deltatheta_gamma = log(mean_deltatheta_gamma); %p-value = 0.029 ONLY WITHOUT CONSEC AMP REQUIREMENT
% 
% coh_delta_theta = coh_delta_theta(~isnan(mean_deltatheta_gamma)); %remove nonsig subj 
% coh_delta_theta = log(coh_delta_theta);
% 
% lm = fitlm(mean_deltatheta_gamma,coh_delta_theta)
% rsq = lm.Rsquared.Adjusted;
% p_val = table2array(lm.Coefficients(2,4));
% r=corr(mean_deltatheta_gamma.',coh_delta_theta.','type','pearson');
% 
% 
% scatter(mean_deltatheta_gamma,coh_delta_theta)
% ylabel('Mean Delta-Theta Coherence')
% xlabel('Mean Delta-Theta Gamma PLV')
% lsline
% annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))
% 

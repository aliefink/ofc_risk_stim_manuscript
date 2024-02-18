
%% Data Loading


%cfc_all_subj = get_cfc_results_all(subjs,data_path); 
%save([data_path 'cfc_results/cfc_all_subj.mat'], 'cfc_all_subj')  
load([data_path 'cfc_results/cfc_all_subj.mat'])
% load([data_path 'cfc_results/cfc_all_subj_sig_info.mat'])
% load([data_path 'cfc_results/sig_elecs_by_subj.mat'])
% load([data_path 'cfc_results/nonsig_elecs_by_subj.mat'])
% load([data_path 'cfc_results/find_sig_elecs_header.mat'])

[sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'prop','no'); % p value calculations from proportion perms > pixel/num perms
[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);
% save([data_path 'cfc_results/sig_mean_plv.mat'], 'sig_mean_plv')%https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf

%brewermap(3,'Set1'); < this is how you get cbrewer palettes!!

%% Extract mean plv data by freq band
[cfc_all_subj] = get_summary_stats(cfc_all_subj, subjs);


[all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
    get_plv_freq_block_method(cfc_all_subj, subjs);

%% Figure settings 

saez_cmap = {'#004v67E', '#561C46', '#FFAC4D', '#C19BD3','#CDB6B6'}; 

%% Fig 2A - Single electrode PLV plot (subj 6 elec 11)
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
'TickLength',[0 0],'FontName','Arial');
% title({'\rmPhase-Amplitude Coupling Within','Single OFC Electrode'},'FontSize',30,'FontName','Arial');
xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Arial');
ylabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Arial');

% cticks = [-1,1,3,5,7,9,11,13];
% clabels = {'-1','1','3','5','7','9','11','13'};
h = colorbar;
% h.Ticks = cticks;
% h.TickLabels = clabels;
set(get(h, 'ylabel'), 'string', '\fontsize{26}PLV Z-Score','FontSize',22,'FontName','Arial');
set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path subj_id '_' elec_id '_single_elec_plv.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

% saveas(plv_single_elec,[fig_path, num2str(elec_num), '_single_elec_plv.pdf'])

%%%%%summary stats for fig a:
% mean(plv_mat,'all')
% std(plv_mat,[],'all','omitnan')
% max(plv_mat,[],'all')
% [row,col] = find(plv_mat == max(plv_mat,[],'all')) %row = amp idx, col = phase idx
% phase_f_array(col)
% amp_f_array(row)

    

%% Fig 2B - Average PLV significant electrodes

sig_plv_group = figure('Name','mean zscore PLV significant elecs across subj')
sig_plv_group.Position = [0, 0, 750, 600]; %best fig size!
contourf(sig_mean_plv);
amp_f_array = 5:5:200;
phase_f_array = 2:20;
% Set locations of ticks and their labels.
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end),...
'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end),'FontSize',22,...
'TickLength',[0 0],'FontName','Arial');
% title({'\rmOFC Phase-Amplitude Coupling','Across Subjects'},'FontSize',30,'FontName','Arial');
% title({'\rmStrong Coupling Between Lower Frequency','Phases & Higher Frequency Amplitudes', 'Consistent Across Subjects'},'FontSize',30,'FontName','Arial');
xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Arial');
ylabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Arial');
h = colorbar;
% h.Ticks = cticks;
% h.TickLabels = clabels;
%change contour visual with clim([-20 20])
set(get(h, 'ylabel'), 'string', '\fontsize{26}PLV Z-Score','FontSize',22,'FontName','Arial');
set(gcf,'PaperOrientation','landscape'); %set(gca, 'box', 'off')


% print(gcf,[fig_path 'group_mean_plv_sig_elecs_short_title.pdf'],'-dpdf','-bestfit')

%%%%%summary stats for fig b:
% mean(sig_mean_plv,'all')
% std(sig_mean_plv,[],'all','omitnan')
% max(sig_mean_plv,[],'all')
% [row,col] = find(sig_mean_plv == max(sig_mean_plv,[],'all')) %row = amp idx, col = phase idx
% phase_f_array(col)
% amp_f_array(row)


% %%%%get contour distribution for max plvs 
% [cm,h] = contourf(sig_mean_plv);
% [contourTable, contourArray] = getContourLineCoordinates(cm);
% max_z = 3;%max(contourTable.Level);
% max_level = contourTable.Level>=max_z;
% max_levelX = contourTable.X(max_level);
% max_levelY = contourTable.Y(max_level); %amplitudes of most significant cluster
% amp_f_array = 5:5:200;
% phase_f_array = 2:20;
% sig_phase_idx = unique(round(max_levelX));
% sig_amp_idx = unique(round(max_levelY));
% sig_phases = phase_f_array(sig_phase_idx');
% sig_amps = amp_f_array(sig_amp_idx');
% sig_phases
% sig_amps
%% Fig 2C - Line plots of frequency distributions
phase_f_array = 2:20;
figure
fig = gcf;

options.handle = fig;
options.color_area = [0 70 126]./255; %6th of 7-class YlGnBu colorbrewer2  33,113,181 33,113,181
options.color_line = [0 70 126]./255;  %8,69,148 8,69,148
options.alpha      = 0.5;
options.line_width = 4;
options.error      = 'sem';
plot_areaerrorbar(sig_mean_plv, options);
hold on 
% set(fig,'defaultAxesFontSize',18);
set(gca,'FontSize',16,'XTick', phase_f_array(1:2:end),'TickLength',[0 0],...
    'FontName','Arial','box','off');
xlabel('Frequency for Phase (Hz)','FontSize',24,'FontName','Arial');
ylabel('PLV Z-Score','FontSize',24,'FontName','Arial');
% ylabel('PLV_{z}','FontSize',24,'FontName','Arial');
% xlim([phase_f_array(1),phase_f_array(end)])
ylim([1 3.5])
hold on

% print(gcf,[fig_path 'GREEN Phase Freq Distribution of PLVs.pdf'],'-dpdf','-bestfit')

plv_distr_phase_freqs = mean(sig_mean_plv,1); %mean along columns - phase means
[phase_peaks, phase_locs] = findpeaks(plv_distr_phase_freqs); % PHASE PEAKS AT 2 AND 4!
phase_peak_yidx = max(plv_distr_phase_freqs) + 0.25;

d = [2 3];
t = [4 5 6 7 8];

for x=1:length(phase_locs)
    freq = round(phase_locs(x));
    if sum(find(d ==freq)) 
    text(freq,phase_peak_yidx,'\fontsize{22}\fontname{Arial}\bf\delta','horizontalalignment','center')
    else
    text(freq,phase_peak_yidx,'\fontsize{22}\fontname{Arial}\bf\theta','horizontalalignment','center')
    end
end 
% print(gcf,[fig_path 'Phase Freq Distribution of PLVs.pdf'],'-dpdf','-bestfit')

% %%%summary stats
% phase_peaks %PLV for each peak
% phase_locs %amp for each peak
% phase_f_array(phase_locs(1))
% phase_f_array(phase_locs(2))
% mean(plv_distr_phase_freqs)
% std(plv_distr_phase_freqs)
% (std(plv_distr_phase_freqs)./sqrt(size(plv_distr_phase_freqs',1))) %SEM


%% AMP

amp_f_array = 5:5:200;
figure
fig = gcf;

options.handle = fig;
options.color_area = [86 28 70]./255; % 239,59,446th of 7-class Reds colorbrewer2 O 217,72,1   R 239,59,44  G 0,109,44 P 106,81,163
options.color_line = [86 28 70]./255; % 153,0,13 166,54,3    153,0,13  166,54,3 0,68,27 84,39,143
options.alpha      = 0.5;  
options.line_width = 4;
options.error      = 'sem'; % 95 percent confidence interval
plot_areaerrorbar(sig_mean_plv', options);
hold on 
% set(fig,'defaultAxesFontSize',18);
set(gca,'XTick',4:4:40,'XTickLabel',amp_f_array(4:4:end),'FontSize',16,...
    'TickLength',[0 0],'FontName','Arial','box','off');
xlabel('Frequency for Amplitude (Hz)','FontSize',24,'FontName','Arial');
% ylabel('PLV_{z}','FontSize',26,'FontName','Arial');
ylabel('PLV Z-Score','FontSize',24,'FontName','Arial');
% xlim([1,40])
hold on 
plv_distr_amp_freqs = mean(sig_mean_plv,2)'; %mean along rows - amp means
[amp_peaks, amp_locs] = findpeaks(plv_distr_amp_freqs); % PHASE PEAKS AT 2 AND 4!


b = 12:1:29;
g = 30:1:59;
hg = 60:1:200;

for x=1:length(amp_locs)
    freq_idx = round(amp_locs(x));
    freq = amp_f_array(round(amp_locs(x)));
    amp_peak_yidx = plv_distr_amp_freqs(1,freq_idx) + 0.5;
    if sum(find(b ==freq)) 
    text(freq_idx,amp_peak_yidx,...
        '\fontsize{20}\fontname{Arial}\bf\beta','horizontalalignment','center')
    elseif sum(find(g ==freq)) 
    text(freq_idx,amp_peak_yidx,...
        '\fontsize{20}\fontname{Arial}\bf\gamma','horizontalalignment','center')
    else
    text(freq_idx,amp_peak_yidx,...
        '\fontsize{20}\fontname{Arial}\itH\rm\fontsize{20}\fontname{Arial}\bf\gamma',...
        'horizontalalignment','center')
    end
    
end


ylim([0,max(plv_distr_amp_freqs +1)]);
% print(gcf,[fig_path 'Amp Freq Distribution of PLVs.pdf'],'-dpdf','-bestfit')

%%%summary stats
amp_peaks %PLV for each peak
amp_locs %amp for each peak
amp_f_array(amp_locs(1))
amp_f_array(amp_locs(2))
mean(plv_distr_amp_freqs)
std(plv_distr_amp_freqs)
(std(plv_distr_amp_freqs)./sqrt(size(plv_distr_amp_freqs',1))) %SEM


%%%%AMP PEAKS AT 35 AND 90

% [pks,locs,w,p] = findpeaks(plv_distr_amp_freqs) 
% findpeaks(plv_distr_amp_freqs)this will find the peaks right on plot!


%to do with just broadgamma - 
% broadgamma = [30 200];
% broadgamma_idx = find((amp_f_array>=broadgamma(1)&amp_f_array<=broadgamma(2)));
% broadgamma_plvs_mean = sig_mean_plv(broadgamma_idx,:);
% 
% broadgamma_distr_phase_freqs = mean(broadgamma_plvs_mean,1);
% 
% plot(phase_f_array,broadgamma_distr_phase_freqs)

%% Join phase/amp dist plots
amp_f_array = 5:5:200;
phase_f_array = 2:20;

freq_dist = figure('Name','Frequency Distributions');
freq_dist.Position = [0,0,700, 600];

tiledlayout(2,1)

% PHASE PLOT
nexttile
options.handle = freq_dist;
options.color_area = [0 70 126]./255; %6th of 7-class YlGnBu colorbrewer2  33,113,181 33,113,181
options.color_line = [0 70 126]./255;  %8,69,148 8,69,148
options.alpha      = 0.5;
options.line_width = 4;
options.error      = 'sem';
plot_areaerrorbar(sig_mean_plv, options);
hold on 
% set(fig,'defaultAxesFontSize',18);
set(gca,'FontSize',16,'XTick', phase_f_array(1:2:end),'TickLength',[0 0],...
    'FontName','Arial','box','off');
xlabel('Frequency for Phase (Hz)','FontSize',24,'FontName','Arial');
ylabel('PLV Z-Score','FontSize',24,'FontName','Arial');
% ylabel('PLV_{z}','FontSize',24,'FontName','Arial');
xlim([1,phase_f_array(end-1)])
% ylim([1 3.5])
hold on

% print(gcf,[fig_path 'GREEN Phase Freq Distribution of PLVs.pdf'],'-dpdf','-bestfit')

plv_distr_phase_freqs = mean(sig_mean_plv,1); %mean along columns - phase means
[phase_peaks, phase_locs] = findpeaks(plv_distr_phase_freqs); % PHASE PEAKS AT 2 AND 4!
phase_peak_yidx = max(plv_distr_phase_freqs) + 0.35;

d = [2 3];
t = [4 5 6 7 8];

for x=1:length(phase_locs)
    freq = round(phase_locs(x));
    if sum(find(d ==freq)) 
    text(freq,phase_peak_yidx,'\fontsize{22}\fontname{Arial}\bf\delta','horizontalalignment','center')
    else
    text(freq,phase_peak_yidx,'\fontsize{22}\fontname{Arial}\bf\theta','horizontalalignment','center')
    end
end 
ylim([0,max(plv_distr_phase_freqs +0.75)]);



% AMP PLOT
nexttile
options.handle = freq_dist;
options.color_area = [86 28 70]./255; % 239,59,446th of 7-class Reds colorbrewer2 O 217,72,1   R 239,59,44  G 0,109,44 P 106,81,163
options.color_line = [86 28 70]./255; % 153,0,13 166,54,3    153,0,13  166,54,3 0,68,27 84,39,143
options.alpha      = 0.5;  
options.line_width = 4;
options.error      = 'sem'; % 95 percent confidence interval
plot_areaerrorbar(sig_mean_plv', options);
hold on 
% set(fig,'defaultAxesFontSize',18);
set(gca,'XTick',4:4:40,'XTickLabel',amp_f_array(4:4:end),'FontSize',16,...
    'TickLength',[0 0],'FontName','Arial','box','off');
xlabel('Frequency for Amplitude (Hz)','FontSize',24,'FontName','Arial');
% ylabel('PLV_{z}','FontSize',26,'FontName','Arial');
ylabel('PLV Z-Score','FontSize',24,'FontName','Arial');
% xlim([1,40])
hold on 
plv_distr_amp_freqs = mean(sig_mean_plv,2)'; %mean along rows - amp means
[amp_peaks, amp_locs] = findpeaks(plv_distr_amp_freqs); % PHASE PEAKS AT 2 AND 4!


b = 12:1:29;
g = 30:1:59;
hg = 60:1:200;

for x=1:length(amp_locs)
    freq_idx = round(amp_locs(x));
    freq = amp_f_array(round(amp_locs(x)));
    amp_peak_yidx = plv_distr_amp_freqs(1,freq_idx);
    if sum(find(b ==freq)) 
    text(freq_idx,amp_peak_yidx,...
        '\fontsize{22}\fontname{Arial}\bf\beta','horizontalalignment','center')
    elseif sum(find(g ==freq)) 
    text(freq_idx,amp_peak_yidx+ 0.65,...
        '\fontsize{22}\fontname{Arial}\bf\gamma','horizontalalignment','center')
    else
    text(freq_idx,amp_peak_yidx+ 0.5,...
        '\fontsize{20}\fontname{Arial}\itH\rm\fontsize{22}\fontname{Arial}\bf\gamma',...
        'horizontalalignment','center')
    end
    
end
ylim([0,max(plv_distr_amp_freqs +1)]);
xlim([1,40])
% set(gcf,'PaperOrientation','landscape'); %set(gca, 'box', 'off')

set(gca, 'box', 'off')



% print(gcf,[fig_path 'Freq Distribution Combo Vert.pdf'],'-dpdf','-bestfit')


%% Risk correlation scatter plot DELTA THETA GAMMA RISK

% mean_deltatheta_gamma = all_elec_plv_info.delta_theta_hga_means;
mean_deltatheta_gamma = sig_elec_plv_info.delta_theta_broadgamma_means(~isnan(sig_elec_plv_info.delta_theta_broadgamma_means)); %SIG
%mean_deltatheta_gamma = log(mean_deltatheta_gamma);
% nonsig_subj_idx = find(isnan(mean_deltatheta_gamma));
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
set(gca, 'box', 'off')
% set(gcf,'PaperOrientation','landscape');
% print(gcf,[fig_path 'plv_risk_broadgamma.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

%without log p-value = 0.0355, with log p-value = 0.147

%delta theta gamma p 0.366, hga p 0.22 broadgamma p 0.249
%alpha p 0.554 beta p 0.419


% outlier_idx = find(mean_deltatheta_gamma == max(mean_deltatheta_gamma));
% mean_deltatheta_gamma(outlier_idx) = [];
% risk_prefs(outlier_idx) = [];



%% Histogram of permutation values 

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

%% FDR Correction Schematic

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
histogram(og_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05);%,'DisplayName','Uncorrected p values');
hold on;
histogram(fdr_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05);%,'DisplayName','FDR Adjusted p values\newline num sig'); %set bin width to alpha for count of < alpha
hold on;
xline(0.05,'--k',...
    '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
    'LineWidth',1.5,'Alpha',1,...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
ylabel('Count','FontSize',22,'FontName','Arial');
% xlim([0,max(fdr_pval_vec)])
legend({join(['Uncorrected \itp\rm values \newline' num_og_sig '/' num_pixels ' < \fontsize{20}\bf\alpha'],''),...
    join(['FDR Corrected \itp\rm values \newline' num_fdr_sig '/' num_pixels ' < \fontsize{20}\bf\alpha'],'')},...
    'FontSize',18,'FontName','Arial','Location','best')
set(gca, 'box', 'off')

set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  '/SOM_methods/significant_pval_two_dist.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

%% Separate the histograms
%uncorrected pvalues
og_pvals = data.prop_pvals{1,find(data.ofc_elecs==elec_num)};
og_pval_vec = reshape(og_pvals,1,[]);
num_og_sig = string(sum(og_pval_vec<0.05));

pval_fig = figure('Name','pvalue distributions') 
pval_fig.Position = [0, 0, 750, 600]; %best fig size!
histogram(og_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05,...
    'DisplayName',join(['\fontname{Arial}\fontsize{18}Uncorrected \itp\rm values \newline '...
    num_og_sig '/' num_pixels ' < \fontsize{22}\bf\alpha\rm\fontsize{18}'],''));
hold on;
xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
    '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
ylabel('Count','FontSize',22,'FontName','Arial');
legend()
set(gca, 'box', 'off')
set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  '/SOM_methods/significant_pval_uncorrected_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page



%corrected pvalues
fdr_pvals = sig_elecs_by_subj{6,find(strcmp(subjs,subj_id))}{find(sig_elecs_by_subj{2,find(strcmp(subjs,subj_id))}==elec_num),1}; %adj pval is row 6,subj is col, elec is row in cell 
fdr_pval_vec = reshape(fdr_pvals,1,[]);
num_fdr_sig = string(sum(fdr_pval_vec<0.05));
%total number of pixels
num_pixels = string(length(fdr_pval_vec));

pval_fig = figure('Name','pvalue distributions') 
pval_fig.Position = [0, 0, 750, 600]; %best fig size!
histogram(fdr_pval_vec,'FaceColor','#D95319','FaceAlpha',0.8,...
    'Normalization','count','BinWidth',0.05,'DisplayName',...
    join(['FDR Corrected \itp\rm values \newline' num_fdr_sig ...
    '/' num_pixels ' < \fontsize{22}\bf\alpha'],'')); %set bin width to alpha for count of < alpha
hold on;
xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
    '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
ylabel('Count','FontSize',22,'FontName','Arial');
% xlim([0,max(fdr_pval_vec)])
legend()

set(gca, 'box', 'off')

set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  '/SOM_methods/significant_pval_fdr_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

%% %% NONSIGFDR Correction Schematic sig = s06 e 11, nonsig = s02 e 22
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
histogram(og_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05);%,'DisplayName','Uncorrected p values');
hold on;
histogram(fdr_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05);%,'DisplayName','FDR Adjusted p values\newline num sig'); %set bin width to alpha for count of < alpha
hold on;
xline(0.05,'--k',...
    '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
    'LineWidth',1.5,'Alpha',1,...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
ylabel('Count','FontSize',22,'FontName','Arial');
% xlim([0,max(fdr_pval_vec)])
legend({join(['Uncorrected \itp\rm values \newline' num_og_sig '/' num_pixels ' < \fontsize{20}\bf\alpha'],''),...
    join(['FDR Corrected \itp\rm values \newline' num_fdr_sig '/' num_pixels ' < \fontsize{20}\bf\alpha'],'')},...
    'FontSize',18,'FontName','Arial','Location','best')
set(gca, 'box', 'off')

set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  '/SOM_methods/nonsignificant_pval_two_dist_s2e21.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page


pval_fig = figure('Name','pvalue distributions') 
pval_fig.Position = [0, 0, 750, 600]; %best fig size!
histogram(og_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05,...
    'DisplayName',join(['\fontname{Arial}\fontsize{18}Uncorrected \itp\rm values \newline '...
    num_og_sig '/' num_pixels ' < \fontsize{22}\bf\alpha\rm\fontsize{18}'],''));
hold on;
xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
    '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
ylabel('Count','FontSize',22,'FontName','Arial');
legend()
set(gca, 'box', 'off')
set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  '/SOM_methods/NONsignificant_pval_uncorrected_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page



pval_fig = figure('Name','pvalue distributions') 
pval_fig.Position = [0, 0, 750, 600]; %best fig size!
histogram(fdr_pval_vec,'FaceColor','#D95319','FaceAlpha',0.8,...
    'Normalization','count','BinWidth',0.05,'DisplayName',...
    join(['FDR Corrected \itp\rm values \newline' num_fdr_sig ...
    '/' num_pixels ' < \fontsize{22}\bf\alpha'],'')); %set bin width to alpha for count of < alpha
hold on;
xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
    '\fontname{Arial}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Arial}\fontsize{18} = 0.05',...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Arial');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Arial');
ylabel('Count','FontSize',22,'FontName','Arial');
% xlim([0,max(fdr_pval_vec)])
legend()

set(gca, 'box', 'off')

set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  '/SOM_methods/NONsignificant_pval_fdr_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page




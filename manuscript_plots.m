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
% load([data_path 'cfc_results/cfc_all_subj_sig_info.mat'])
% load([data_path 'cfc_results/sig_elecs_by_subj.mat'])
% load([data_path 'cfc_results/nonsig_elecs_by_subj.mat'])
% load([data_path 'cfc_results/find_sig_elecs_header.mat'])

[sig_elecs_by_subj,nonsig_elecs_by_subj,sig_count,nonsig_count, find_sig_elecs_header,cfc_all_subj] =find_sig_elecs(cfc_all_subj,subjs,'prop','no'); % p value calculations from proportion perms > pixel/num perms
[sig_mean_plv,all_sig_plvs,nonsig_mean_plv,all_nonsig_plvs] = plot_mean_phase_amp_ranges(cfc_all_subj, subjs, sig_count,nonsig_count);
% save([data_path 'cfc_results/sig_mean_plv.mat'], 'sig_mean_plv')%https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf

%brewermap(3,'Set1'); < this is how you get cbrewer palettes!!

%% Extract mean plv data by freq band


[all_elec_plv_info,sig_elec_plv_info,nonsig_elec_plv_info,phase_freqs,amp_freqs] = ...
    get_plv_freq_block_method(cfc_all_subj, subjs);

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
'TickLength',[0 0],'FontName','Lucida Grande');
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
'TickLength',[0 0],'FontName','Lucida Grande');
title({'\rmOFC Phase-Amplitude Coupling','Across Subjects'},'FontSize',30,'FontName','Lucida Grande');
% title({'\rmStrong Coupling Between Lower Frequency','Phases & Higher Frequency Amplitudes', 'Consistent Across Subjects'},'FontSize',30,'FontName','Lucida Grande');
xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Lucida Grande');
ylabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Lucida Grande');
h = colorbar;
% h.Ticks = cticks;
% h.TickLabels = clabels;
%change contour visual with clim([-20 20])
set(get(h, 'ylabel'), 'string', '\fontsize{26}PLV Z-Score','FontSize',22,'FontName','Lucida Grande');
set(gcf,'PaperOrientation','landscape','box','off'); %set(gca, 'box', 'off')


% print(gcf,[fig_path 'group_mean_plv_sig_elecs_short_title.pdf'],'-dpdf','-bestfit')

%% Fig 2C - Line plots of frequency distributions
phase_f_array = 2:20;
figure
fig = gcf;

options.handle = fig;
options.color_area = [33,113,181]./255; %6th of 7-class YlGnBu colorbrewer2  33,113,181
options.color_line = [8,69,148]./255;  8,69,148
options.alpha      = 0.5;
options.line_width = 4;
options.error      = 'sem';
plot_areaerrorbar(sig_mean_plv, options);
hold on 
% set(fig,'defaultAxesFontSize',18);
set(gca,'FontSize',16,'XTick', phase_f_array(1:2:end),'TickLength',[0 0],...
    'FontName','Helvetica','box','off');
xlabel('Frequency for Phase (Hz)','FontSize',24,'FontName','Helvetica');
ylabel('PLV Z-Score','FontSize',24,'FontName','Helvetica');
% ylabel('PLV_{z}','FontSize',24,'FontName','Lucida Grande');
% xlim([phase_f_array(1),phase_f_array(end)])
% ylim([1 3.25])
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
    text(freq,phase_peak_yidx,'\fontsize{20}\fontname{Helvetica}\bf\delta\cdot\gamma','horizontalalignment','center')
    else
    text(freq,phase_peak_yidx,'\fontsize{20}\fontname{Helvetica}\bf\theta\cdot\gamma','horizontalalignment','center')
    end
end 


%%%%%PHASE line plot

% % phase_f_array = 2:20;
% %phase frequency distribution by mean plvs
% figure('Name','Phase freq PLV Dist')
% plv_distr_phase_freqs = mean(sig_mean_plv,1); %mean along columns - phase means
% plot(phase_f_array,plv_distr_phase_freqs,'LineWidth',4,'Color','#225ea8') %2171b5  225ea8
% set(gca,'FontSize',18,'XTick', phase_f_array(1:2:end),'TickDir','in','FontName','Lucida Grande','box','off');
% xlabel('Frequency for Phase (Hz)','FontSize',26,'FontName','Lucida Grande');
% ylabel('PLV_{z}','FontSize',26,'FontName','Lucida Grande');
% xlim([phase_f_array(1),phase_f_array(end)])
% ylim([1 3.25])
% %add error bars 
% hold on
% phase_sem = std(sig_mean_plv)/sqrt(19);
% errorbar(phase_f_array,plv_distr_phase_freqs,phase_sem,'k','linestyle','none','LineWidth', 1.5,'CapSize',1,'Color','#225ea8');%'Color', [.25 .25 .25]);
%% AMP

amp_f_array = 5:5:200;
figure
fig = gcf;

options.handle = fig;
options.color_area = [239,59,44]./255; %6th of 7-class Reds colorbrewer2 O 217,72,1   R 239,59,44  G 0,109,44 P 106,81,163
options.color_line = [153,0,13]./255; %166,54,3    153,0,13  166,54,3 0,68,27 84,39,143
options.alpha      = 0.5;  
options.line_width = 4;
options.error      = 'sem'; % 95 percent confidence interval
plot_areaerrorbar(sig_mean_plv', options);
hold on 
% set(fig,'defaultAxesFontSize',18);
set(gca,'XTick',4:4:40,'XTickLabel',amp_f_array(4:4:end),'FontSize',16,...
    'TickLength',[0 0],'FontName','Lucida Grande','box','off');
xlabel('Frequency for Amplitude (Hz)','FontSize',24,'FontName','Lucida Grande');
% ylabel('PLV_{z}','FontSize',26,'FontName','Lucida Grande');
ylabel('PLV Z-Score','FontSize',24,'FontName','Lucida Grande');
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
        '\fontsize{20}\fontname{Helvetica}\bf\beta','horizontalalignment','center')
    elseif sum(find(g ==freq)) 
    text(freq_idx,amp_peak_yidx,...
        '\fontsize{20}\fontname{Helvetica}\bf\gamma','horizontalalignment','center')
    else
    text(freq_idx,amp_peak_yidx,...
        '\fontsize{20}\fontname{Helvetica}\itH\rm\fontsize{20}\fontname{Helvetica}\bf\gamma',...
        'horizontalalignment','center')
    end
    
end


ylim([0,max(plv_distr_amp_freqs +1)]);

% print(gcf,[fig_path 'Amp Freq Distribution of PLVs.pdf'],'-dpdf','-bestfit')
%%%%AMP PEAKS AT 35 AND 90

% %% AMP line plot
% amp_f_array = 5:5:200;
% %amp frequency distribution by mean plvs
% figure('Name','Amp freq PLV Dist')
% plv_distr_amp_freqs = mean(sig_mean_plv,2)'; %mean along rows - amp means
% plot(amp_f_array,plv_distr_amp_freqs,'LineWidth',4,'Color','#cb181d') %max peak @ 35 second @ 90
% set(gca,'XTick',amp_f_array(4:4:end),'FontSize',18,'TickDir','in','FontName','Lucida Grande','box','off');
% xlabel('Frequency for Amplitude (Hz)','FontSize',26,'FontName','Lucida Grande');
% ylabel('PLV_{z}','FontSize',26,'FontName','Lucida Grande');
% xlim([amp_f_array(1),amp_f_array(end)])
% ylim([0 3.25])
% %add error bars 
% %calculate SEM for each amp freq
% hold on 
% amp_sem = std(sig_mean_plv')/sqrt(40);
% errorbar(amp_f_array,plv_distr_amp_freqs,amp_sem,'k','linestyle','none','LineWidth', 1.5,'CapSize',1,'Color','#cb181d');%'Color', [.25 .25 .25]);

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


%% SIG ELEC ONLY BARPLOT - BROADGAMMA
%%%%%   BROADGAMMA ONLY

plv_by_freq = figure('Name','Mean PLV by Freq');
plv_by_freq.Position = [0,0,700, 600];% plv_by_freq.Position(3:4) = [1000 1000];

x_idxs = [1 2 3];
bar_y_vals = [log(mean(sig_elec_plv_info.delta_theta_broadgamma_means,'omitnan'));...
    log(mean(sig_elec_plv_info.alpha_broadgamma_means,'omitnan'));...
    log(mean(sig_elec_plv_info.beta_broadgamma_means,'omitnan'))];
% 
% bar_y_vals = [mean(sig_elec_plv_info.delta_theta_broadgamma_means,'omitnan');...
%     mean(sig_elec_plv_info.alpha_broadgamma_means,'omitnan');...
%     mean(sig_elec_plv_info.beta_broadgamma_means,'omitnan')];

b = bar(x_idxs,bar_y_vals,'FaceColor','flat','EdgeColor', 'w');

%words only XTickLabels

% 


set(gca,'FontSize',18,'FontName','Lucida Grande','TickLength',[0 0])
xticklabels({'\fontsize{28}\fontname{Lucida Grande}\bf\delta\theta\cdot\gamma',...
    '\fontsize{28}\fontname{Lucida Grande}\bf\alpha\cdot\gamma',...
    '\fontsize{28}\fontname{Lucida Grande}\bf\beta\cdot\gamma'})
% xticklabels({'\fontsize{24}Delta-Theta',...
%     '\fontsize{24}Alpha',...
%     '\fontsize{24}Beta'})
xlabel('\fontname{Lucida Grande}Phase-Amplitude Frequencies (Hz)','FontSize',26)
ylabel('\fontname{Lucida Grande}PLV Z-Score','FontSize',26)


% title({'\rm\fontname{Lucida Grande}Sig Elec Coupling Between Low Frequency Phase and'; 'Broadband Gamma Amplitude'},...
%     'VerticalAlignment','bottom','FontSize',24)


hold on 


%manually set color data for each bar
b(1,1).CData(1,:) = [203,24,29]./255;%[203,24,29]./255;brightR[203,24,29]deepR[165,15,21]Or[166,54,3]
b(1,1).CData(2,:) = [252,146,114]./255;%[252,146,114]./255;Or[241,105,19]
b(1,1).CData(3,:) = [254,224,210]./255;%[254,224,210]./255;Or[253,174,107]

hold on 
%add sem bars to amp frequency label to top of bars
%calculate sem for bar
b1_data_pts = log(sig_elec_plv_info.delta_theta_broadgamma_means);
b1_sem =  log(std(sig_elec_plv_info.delta_theta_broadgamma_means,'omitnan'))/sqrt(14);
b1_upper = log(mean(b1_data_pts,'omitnan'))+b1_sem;
b1_lower = log(mean(b1_data_pts,'omitnan'))-b1_sem;

b2_data_pts = log(sig_elec_plv_info.alpha_broadgamma_means);
b2_sem = log(std(sig_elec_plv_info.alpha_broadgamma_means,'omitnan'))/sqrt(14);
b2_upper = log(mean(b2_data_pts,'omitnan'))+b2_sem;
b2_lower = log(mean(b2_data_pts,'omitnan'))-b2_sem;

b3_data_pts = log(sig_elec_plv_info.beta_broadgamma_means);
b3_sem =  log(std(sig_elec_plv_info.beta_broadgamma_means,'omitnan'))/sqrt(14);
b3_upper = log(mean(b3_data_pts,'omitnan'))+b3_sem;
b3_lower = log(mean(b3_data_pts,'omitnan'))-b3_sem;
% 
% b1_data_pts = log(sig_elec_plv_info.delta_theta_broadgamma_means);
% b1_sem =  log(std(sig_elec_plv_info.delta_theta_broadgamma_means,'omitnan'))/sqrt(14);
% b1_upper = log(mean(b1_data_pts,'omitnan'))+b1_sem;
% b1_lower = log(mean(b1_data_pts,'omitnan'))-b1_sem;
% 
% b2_data_pts = log(sig_elec_plv_info.alpha_broadgamma_means);
% b2_sem = log(std(sig_elec_plv_info.alpha_broadgamma_means,'omitnan'))/sqrt(14);
% b2_upper = log(mean(b2_data_pts,'omitnan'))+b2_sem;
% b2_lower = log(mean(b2_data_pts,'omitnan'))-b2_sem;
% 
% b3_data_pts = log(sig_elec_plv_info.beta_broadgamma_means);
% b3_sem =  log(std(sig_elec_plv_info.beta_broadgamma_means,'omitnan'))/sqrt(14);
% b3_upper = log(mean(b3_data_pts,'omitnan'))+b3_sem;
% b3_lower = log(mean(b3_data_pts,'omitnan'))-b3_sem;

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
errorbar(xtips1',err_data,err_sem,'k','linestyle','none','LineWidth', 1.5,'CapSize',5);%'Color', [.25 .25 .25]);


set(gca, 'box', 'off')
set(gcf,'PaperOrientation','landscape');
% print(gcf,[fig_path 'greeksymbols barplot_broadgamma_sem.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page

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


%significance testing for barplot 


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
    'TickLength',[0 0],'FontSize',22,'FontName','Lucida Grande')%...

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
% set(gca, 'box', 'off')






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

%% FDR Correction Schematic sig = s06 e 11, nonsig = s02 e 22
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
histogram(og_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05);%,'DisplayName','Uncorrected p values');
hold on;
histogram(fdr_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05);%,'DisplayName','FDR Adjusted p values\newline num sig'); %set bin width to alpha for count of < alpha
hold on;
xline(0.05,'--k',...
    '\fontname{Lucida Grande}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Lucida Grande}\fontsize{18} = 0.05',...
    'LineWidth',1.5,'Alpha',1,...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Lucida Grande');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Lucida Grande');
ylabel('Count','FontSize',22,'FontName','Lucida Grande');
% xlim([0,max(fdr_pval_vec)])
legend({join(['Uncorrected \itp\rm values \newline' num_og_sig '/' num_pixels ' < \fontsize{20}\bf\alpha'],''),...
    join(['FDR Corrected \itp\rm values \newline' num_fdr_sig '/' num_pixels ' < \fontsize{20}\bf\alpha'],'')},...
    'FontSize',18,'FontName','Lucida Grande','Location','best')
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
    'DisplayName',join(['\fontname{Lucida Grande}\fontsize{18}Uncorrected \itp\rm values \newline '...
    num_og_sig '/' num_pixels ' < \fontsize{22}\bf\alpha\rm\fontsize{18}'],''));
hold on;
xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
    '\fontname{Lucida Grande}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Lucida Grande}\fontsize{18} = 0.05',...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Lucida Grande');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Lucida Grande');
ylabel('Count','FontSize',22,'FontName','Lucida Grande');
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
    '\fontname{Lucida Grande}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Lucida Grande}\fontsize{18} = 0.05',...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Lucida Grande');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Lucida Grande');
ylabel('Count','FontSize',22,'FontName','Lucida Grande');
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
    '\fontname{Lucida Grande}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Lucida Grande}\fontsize{18} = 0.05',...
    'LineWidth',1.5,'Alpha',1,...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Lucida Grande');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Lucida Grande');
ylabel('Count','FontSize',22,'FontName','Lucida Grande');
% xlim([0,max(fdr_pval_vec)])
legend({join(['Uncorrected \itp\rm values \newline' num_og_sig '/' num_pixels ' < \fontsize{20}\bf\alpha'],''),...
    join(['FDR Corrected \itp\rm values \newline' num_fdr_sig '/' num_pixels ' < \fontsize{20}\bf\alpha'],'')},...
    'FontSize',18,'FontName','Lucida Grande','Location','best')
set(gca, 'box', 'off')

set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  '/SOM_methods/nonsignificant_pval_two_dist_s2e21.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page


pval_fig = figure('Name','pvalue distributions') 
pval_fig.Position = [0, 0, 750, 600]; %best fig size!
histogram(og_pval_vec,'FaceAlpha',0.8,'Normalization','count','BinWidth',0.05,...
    'DisplayName',join(['\fontname{Lucida Grande}\fontsize{18}Uncorrected \itp\rm values \newline '...
    num_og_sig '/' num_pixels ' < \fontsize{22}\bf\alpha\rm\fontsize{18}'],''));
hold on;
xline(0.05,'--k','LineWidth',1.5,'Alpha',1,'DisplayName',...
    '\fontname{Lucida Grande}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Lucida Grande}\fontsize{18} = 0.05',...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Lucida Grande');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Lucida Grande');
ylabel('Count','FontSize',22,'FontName','Lucida Grande');
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
    '\fontname{Lucida Grande}\fontsize{18}Significance Threshold\newline\fontsize{22}\bf\alpha\rm\fontname{Lucida Grande}\fontsize{18} = 0.05',...
    'LabelOrientation','horizontal')
set(gca,'TickLength',[0 0],'FontSize',18,'FontName','Lucida Grande');
xlabel('PLV \itp\rm value','FontSize',22,'FontName','Lucida Grande');
ylabel('Count','FontSize',22,'FontName','Lucida Grande');
% xlim([0,max(fdr_pval_vec)])
legend()

set(gca, 'box', 'off')

set(gcf,'PaperOrientation','landscape'); %rotate plot to save horizontally
% print(gcf,[fig_path  '/SOM_methods/NONsignificant_pval_fdr_only.pdf'],'-dpdf','-bestfit') %use print function to save as pdf and 'best fit' the page



%% Risk correlation scatter plot DELTA THETA GAMMA RISK

% mean_deltatheta_gamma = all_elec_plv_info.delta_theta_hga_means;
mean_deltatheta_gamma = sig_elec_plv_info.delta_theta_broadgamma_means(~isnan(sig_elec_plv_info.delta_theta_broadgamma_means)); %SIG
mean_deltatheta_gamma = log(mean_deltatheta_gamma);
% nonsig_subj_idx = find(isnan(mean_deltatheta_gamma));
risk_prefs = [0.4644 0.5474 0.4435 0.529 0.5173 0.539 0.4229 0.3819 0.4914 0.4957 0.5712 0.3628 0.4865 0.4609 0.4853];
risk_prefs = risk_prefs(~isnan(mean_deltatheta_gamma)); %remove nonsig subj 
% risk_prefs = log(risk_prefs);

lm = fitlm(mean_deltatheta_gamma,risk_prefs)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(mean_deltatheta_gamma.',risk_prefs.','type','pearson');

scatter(mean_deltatheta_gamma,risk_prefs)
xlabel('log mean plv btwn delta-theta and hga')
ylabel('Risk Preference')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))

%without log p-value = 0.0355, with log p-value = 0.147

%delta theta gamma p 0.366, hga p 0.22 broadgamma p 0.249
%alpha p 0.554 beta p 0.419
% 

%% Fig 3E - Coherence correlation scatter plot
coh_header = {'RiskPref', 'coh_delta_theta', 'coh_delta',	'coh_theta',	'coh_alpha',  'coh_beta'};
load('/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/data/coh_data/OrG-OrG-Decision_coh_all_lowfreq.mat')
coh_delta_theta = mean_coh(:,3)';
% mean_delta_theta_gamma_plvs = all_elec_plv_info.delta_theta_hga_means; %0.363
mean_deltatheta_gamma = sig_elec_plv_info.delta_theta_broadgamma_means(~isnan(sig_elec_plv_info.delta_theta_broadgamma_means)); %SIG
mean_deltatheta_gamma = log(mean_deltatheta_gamma); %p-value = 0.029 ONLY WITHOUT CONSEC AMP REQUIREMENT

coh_delta_theta = coh_delta_theta(~isnan(mean_deltatheta_gamma)); %remove nonsig subj 
% coh_delta_theta = log(coh_delta_theta);

lm = fitlm(mean_deltatheta_gamma,coh_delta_theta)
rsq = lm.Rsquared.Adjusted;
p_val = table2array(lm.Coefficients(2,4));
r=corr(mean_deltatheta_gamma.',coh_delta_theta.','type','pearson');


scatter(mean_deltatheta_gamma,coh_delta_theta)
ylabel('Mean Delta-Theta Coherence')
xlabel('Mean Delta-Theta Gamma PLV')
lsline
annotation('textbox',[.7 .75 .2 .1],'string',strcat('Adjusted Rsq: ',string(rsq), ' Pearson R: ',string(r)))


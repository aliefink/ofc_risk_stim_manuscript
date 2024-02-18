% Plots electrodes from leGui registration on an MNI template
% Electrode labels are visually curated but follow YBA labels
% Highlights ROIs from electrodes using YBA meshes
% Uses YBA color scheme
% The script can be easily modified to run recon on a single subject
% v1.0 I.S. Apr 2023

%% Prepare
% Load dependencies 
% In addition to these, the script requires elec.csv files and meshes for
% YBA ROIs (pregenerated in Fieldtrip and included in dependency folder)
% export_fig also requires Ghostscript to be installed - get from
% https://pages.uoregon.edu/koch/ (v.10.00 tested and working Apr 2023)

% Requires a couple of dependent functions and Fieldtrip installed
addpath('~cbrewer/')
addpath('~export_fig/')

% My paths

% Load and initialize Fieldtrip
addpath(path_to_fieldtrip)
load([path_to_fieldtrip,'~surface_pial_both.mat']); pial = mesh;
ft_defaults

% Read in ROI table
% This is initially from YBA but modified to include useful ROI definitions
roi = readtable([path_to_script,'/dependencies/yba_roi.csv']);

% Specify your ROI list
roi_list = {'Orbital frontal'}; % ROI list

%% Loop over patients and concatenate in a single elec_all struct
% This requires all SUBID_labels.csv files generated during anatomical recon (LeGui)
% Only electrodes in grey matter are kept
% This is only necessary once - you can save/load elec_all if you have a
% full dataset
label_files = dir([path_to_recons,'*_atlas_loc.csv']);


% Load elec csv, get OFC electrode identity from OrG_elec_names.mat

for i = 1:size(label_files,1)

    % Load electrode csv
    elec = readtable([path_to_recons,label_files(i).name]);
    elec.SUBID = repmat(label_files(i).name(1:5),size(elec,1),1);
    display(['Reading ',label_files(i).name]);

    % Does this subject have OFC elecs?
    %subid = strcat('s',compose("%02.f",i));
    subid = elec.SUBID(1,1:3);
    if(ismember(subid, elec_OrG(:,2))) % If the current subject exists in elec_OrG
        ofc_elecs_temp = elec_OrG(find(contains(elec_OrG(:,2),subid)),1); % Find all elecs in elec_OrG for this subject
        display(['Found ',num2str(size(ofc_elecs_temp,1)), ' elecs from ',subid]);        
        ofc_coord_temp = cat(2, ...
            elec.Coord_x(matches(elec.Electrode,ofc_elecs_temp)), ...
            elec.Coord_y(matches(elec.Electrode,ofc_elecs_temp)), ...
            elec.Coord_z(matches(elec.Electrode,ofc_elecs_temp)));
    end
    display(['Adding ',num2str(size(ofc_coord_temp,1)), ' elecs from ',subid]);

    % Concatenate
    if i == 1; ofc_coord = ofc_coord_temp; else; ofc_coord = [ofc_coord; ofc_coord_temp]; end

end

%% Plot elecs from specific ROIs only
% ROI specification
% Requires reading a yba_roi lookup table
% ROI table has increasingly detailed levels: Hemisphere, Lobe, Region, Gyrus, Code
% ROIs don't need to be at the same level (i.e. can specify Temporal +  Hippocampus)
% In addition there is a Custom colum for manually defined sets of regions (e.g. vmPFC) - relevant ROIs can be defined manually here

% Plot
roi_level = {'Lobe','Region','Gyrus','Code','Custom'};
roi_list = {'Orbital frontal'}; % Example ROI list

fig = figure();
ft_plot_mesh(pial,'facecolor','w','edgecolor',[.5 .5 .5],'facealpha',0.01,'edgealpha',0.05) % a nice looking pial 3d plot
hold on
% Iterate over all ROIs in dpecified list
for i = 1:size(roi_list,2)
    % Find the YBA labels for all regions contained in the current ROI
    roi_level_idx = roi_level{find([contains(roi_list{i},roi.Lobe),contains(roi_list{i},roi.Region),contains(roi_list{i},roi.Gyrus),contains(roi_list{i},roi.Code),contains(roi_list{i},roi.Custom)])}; % Find the right level (Hemisphere/Lobe/Region/Gyrus/Code) for the specified ROI
    roi_code_idx = find(contains(roi.(roi_level_idx),roi_list{i}));
    % For each YBA region, plot the mesh using YBA colors and 3dplot all
    % electrodes contained
    for j = 1:size(roi_code_idx,1)
        idx = roi_code_idx(j);
        load([path_to_script,'/dependencies/YBA/volume_',roi.Code{idx},'_mesh.mat']);
        ft_plot_mesh(mesh,'facecolor','w','edgecolor',[roi.R(idx,:),roi.G(idx,:),roi.B(idx,:)],'facealpha',0.01,'edgealpha',1)
        % Find electrodes in ROI and plot
        %elec_temp = elec_all(find(contains(elec_all.YBA_1,roi.Long_name{idx})),:);
        elec_temp = elec;
    end
end
scatter3(ofc_coord(:,1),ofc_coord(:,2),ofc_coord(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]) % White markers; alternatively use [roi.R(idx,:),roi.G(idx,:),roi.B(idx,:)]) for markers same color as ROI but these are not very visible


fig = figure();
ft_plot_mesh(pial,'facecolor','w','edgecolor',[.5 .5 .5],'facealpha',0.01,'edgealpha',0.05) % a nice looking pial 3d plot
hold on
scatter3(ofc_coord(:,1),ofc_coord(:,2),ofc_coord(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1])


%% Export to pdf (slow!)
set(fig, 'Position',  [100, 100, 1000, 800]) % Sets consistent figure size
view([0 90]); % Sets rotation - ventral view
export_fig(fig,[path_to_recons,'/../../Figures/ofc_ventral.pdf'],'-transparent')
view([-90 0]); % for lateral view
export_fig(fig,[path_to_recons,'/../../Figures/ofc_dorsal.pdf'],'-transparent')

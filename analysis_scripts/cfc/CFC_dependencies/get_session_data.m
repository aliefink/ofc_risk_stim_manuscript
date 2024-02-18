function [chunked_low chunked_high] = get_session_data(SUBID, lf_phase, hf_phase, ...
                                                        cond_string, start_window,...
                                                        end_window)

% Pass subj_globals so buttonpress_times, game_times, etc. can be interpreted
% Specify path - local or cluster
currpath = pwd;
if strcmp(currpath(1:10),'/Users/ign')
    se_path = '/Users/ignaciosaez/Dropbox/Berkeley/Project_EcoG';
elseif strcmp(currpath(1:10),'/home/knig') % Cluster
    se_path = '/home/knight/isaez/Desktop/Project_EcoG';
end
load([se_path '/data/',SUBID,'/subj_globals.mat']);

if ischar(cond_string)
    condition = eval(cond_string);%
elseif isnumeric(cond_string)
    condition = cond_string;
end
%condition = stim_long_only; %for long_trials analysis

% Find time points to keep.
keepers = [];
for i=1:length(condition)
    event_start = round(condition(i)); % Make the time an integer.
    keepers = [keepers (event_start+start_window):(event_start+end_window)];
end

% Make sure there are no repeated time points in keepers (due
% to overlaps of windows with preceding or succeeding trials).
keepers = unique(keepers);

% Make sure there are no indices in keepers beyond the length
% of signal.
while keepers(end) > length(lf_phase)
    keepers(end) = [];
end

chunked_low  = lf_phase(:, keepers);

if ~isempty(hf_phase) %in case get_session_data is only analyzing lf_phase
    chunked_high = hf_phase(:, keepers);
elseif isempty(hf_phase)
    chunked_high = [];
else
    error('chunked_high must either be [] or must exist!')
end

end
function [chunked_low chunked_high] = get_session_data_AF(SUBID, data_path, lf_phase, hf_phase, ...
                                                        cond_string, start_window,...
                                                        end_window)

% AF Edit 9/20/2023 remove trials where button press did not occur 
% Pass subj_globals so buttonpress_times, game_times, etc. can be interpreted
% Specify path - local or cluster
%currpath = pwd;
% if strcmp(currpath(1:10),'/Users/ign')
%     se_path = '/Users/ignaciosaez/Dropbox/Berkeley/Project_EcoG';
% elseif strcmp(currpath(1:10),'/home/knig') % Cluster
%     se_path = '/home/knight/isaez/Desktop/Project_EcoG';
% end
%se_path = '/Users/alexandrafink/Documents/GraduateSchool/SaezLab/gambling_stim_cfc/data/cfc_data/'
load([data_path,'/',SUBID,'_subj_globals.mat']);

bad = find(bad_trials==1); %remove bad trials (time out trials)

if ischar(cond_string)
    condition = eval(cond_string);%
elseif isnumeric(cond_string)
    condition = cond_string;
end
condition(bad)=0; %make bad trials zero index - to remove in next line
condition = nonzeros(condition);
%this is necessary to remove fail trials so function doesn't break


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

chunked_low  = lf_phase(:,keepers); %this line will fail if there are fail trials in time vector

if ~isempty(hf_phase) %in case get_session_data is only analyzing lf_phase
    chunked_high = hf_phase(:, keepers);
elseif isempty(hf_phase)
    chunked_high = [];
else
    error('chunked_high must either be [] or must exist!')
end

end
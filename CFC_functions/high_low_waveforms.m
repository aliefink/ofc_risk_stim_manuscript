function high_low_waveforms(e_mat,srate,subject,cond_string, ...
                            start_window,end_window)

cd(['~/data/' subject '/gdat_CAR_by_elec/'])
load(e_mat)                             % Load electrode_signal.

eeglab
delta_signal = eegfilt(electrode_signal,srate,1,4);
hg_signal = eegfilt(electrode_signal,srate,50,150);

cd('..');
load(cond_string)
condition = round(eval(cond_string));   % Round to make each time an integer.
num_trials = length(condition);
window = start_window:end_window;

% Get delta waveform for average trial.
average_delta = zeros(size(window));
for i=1:num_trials
    event_start = condition(i);
    average_delta = average_delta + delta_signal(condition(i)+window);
end
average_delta = average_delta/num_trials;

% Get hg waveform for average trial.
average_hg = zeros(size(window));
for i=1:num_trials
    event_start = condition(i);
    average_hg = average_hg + hg_signal(condition(i)+window);
end
average_hg = average_hg/num_trials;

close all
figure;
h=subplot(2,1,1);
plot(window, average_delta);
xlim(window([1 end]));
title('1-4 Hz');
h=subplot(2,1,2);
plot(window, average_hg);
xlim(window([1 end]));
title('50-150 Hz');

print('-dpdf', ['~/waveform_figures/' subject '/' cond_string '/' e_mat]);

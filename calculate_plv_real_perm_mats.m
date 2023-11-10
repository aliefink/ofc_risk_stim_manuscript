function [plv_matrix_real, plv_matrix_perm] = calculate_plv_real_perm_mats(data_path,...
    subj_id, cond_string, elec_num, srate, start_window, end_window, num_perm)

tic

%% Initialize

%num_perm  = 2000; % or more; more is better if you have time to run it
%num_perm  = 200; % or more; more is better if you have time to run it

% load signal.mat file and subj_globals.mat
pth_data = [data_path subj_id]; 
load(strcat(pth_data,'/',subj_id,'_signal.mat')) %s13 signal_ds
load(strcat(pth_data,'/',subj_id,'_subj_globals.mat'))


e_signal = signal(elec_num,:);


amp_f_array = 5:5:200;
phase_f_array = 2:20;

% get_signal_parameters, which returns the structure 'sp':
sp = get_signal_parameters('sampling_rate',srate,... % Hz % added srate to function input
    'number_points_time_domain',length(e_signal));

n_amp_freqs = length(amp_f_array); %number of frequencies to analyze for amp data
n_phase_freqs = length(phase_f_array);%numer of frequencies to analyze for phase data
plv_matrix_real = nan(n_amp_freqs, n_phase_freqs);%initialize a matrix of plv values
plv_matrix_perm = nan(n_amp_freqs,n_phase_freqs,num_perm);%initialize a matrix of plv values- this matrix is as deep as the number of permogates run.
% pvalues_matrix = zeros(n_amp_freqs, n_phase_freqs);%initialize a matrix of p-values for each frequency pair

elec = num2str(elec_num);


%% Analysis
% loops over phase (low) frequencies, then by amp (high) frequencies
% For each combination, it calculates and stores a real PLV value and a
% permogate PLV value (after randomly circshifting the phase vector)
for i_phase = 1:n_phase_freqs %for each frequency of the phase data
    
    phase_f = phase_f_array(i_phase); %the particular frequency to filter lf signal
    
    g.center_frequency = phase_f;
    g.fractional_bandwidth = 0.25;
    g.chirp_rate = 0;
    g1 = make_chirplet('chirplet_structure', g, 'signal_parameters', sp);
    
    % filter raw signal at low frequency, extract phase:
    fs = filter_with_chirplet('raw_signal', e_signal, ...
        'signal_parameters', sp, ...
        'chirplet', g1);
    lf_phase = angle(fs.time_domain); %low frequency phase time series %%%should save this out if possible
    clear g.center_frequency
    
    for i_amp = 1:n_amp_freqs %for each frequency of the amplitude data
        
        amp_f = amp_f_array(i_amp); %the particular frequency to filter hf data
        
        g.center_frequency = amp_f;
        g2 = make_chirplet('chirplet_structure', g, 'signal_parameters', sp);
        c
        % filter raw signal at high frequency, extract amplitude:
        fs = filter_with_chirplet('raw_signal', e_signal, ...
            'signal_parameters', sp, ...
            'chirplet', g2);
        hf_amp = abs(fs.time_domain); %%% should save this out at some point!
        
        % filter high frequency amplitude time-series at low frequency, extract phase:
        fs = filter_with_chirplet('raw_signal', hf_amp, ...
            'signal_parameters', sp, ...
            'chirplet', g1);%filter at low frequency
        
        hf_phase = angle(fs.time_domain); %extract phase of high frequency amplitude %%%should save this out!t
        %should plot lf_phase, hf_amp, hf_phase out
        
        % this chunks out data for ALL EVENTS in a particular condition AFTER*** filtering
        [chunked_low_phase chunked_high_phase] = get_session_data_AF(subj_id, pth_data,...
            lf_phase, hf_phase, cond_string, start_window, end_window);
        
        % Compute cross-frequency phase locking value (PLV).
        %plv = abs(mean(exp(1i*(hf_phase - lf_phase))));
        plv = abs(mean(exp(1i*(chunked_high_phase - chunked_low_phase))));
        
        % Store PLV.
        plv_matrix_real(i_amp, i_phase) = plv;
        
        %create permogate distribution of PLV values for each pair of frequencies in the comodulogram
            
        % compute an ensemble of permogate PLVs to compare to actual value to
        % establish statistical significance:
        iteration = [num2str(i_amp) '_' num2str(i_phase)];
        disp(iteration);

        for s = 1:num_perm

            shift = round(rand*sp.number_points_time_domain); %choose a random timepoint
            permogate_lf_phase = circshift(chunked_low_phase,[0 shift]); %shift the timecourse of the low frequency phase data at the random cut point
            plv_matrix_perm(i_amp, i_phase, s) = abs(mean(exp(1i*(chunked_high_phase - permogate_lf_phase))));%permogate plv value
        end            
            
        
    end % end hf amp for loop
    
end % end lf phase for loop
if ~exist([data_path '/cfc_results/' subj_id '/'],'dir') % checks if the appropriate subfolder for this study/block has been created in 'analysis'....
    mkdir([data_path '/cfc_results/' subj_id '/']);
end

% save PLV matrix 
save([data_path '/cfc_results/' subj_id '/' subj_id '_plv_matrix_real_' elec], 'plv_matrix_real');
save([data_path '/cfc_results/' subj_id '/' subj_id '_plv_matrix_perm_' elec], 'plv_matrix_perm');

toc

end

 
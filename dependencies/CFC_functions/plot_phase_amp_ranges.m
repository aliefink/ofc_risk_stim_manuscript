function [plv_matrix] = plot_phase_amp_ranges(SUBID, cond_string, elec_num, ...
    start_window, ...
    end_window, ...
    run_surrogate, ...
    varargin)
tic

% Phase-amplitude cross-frequency coupling analysis
% Creates a comodulogram that plots out the Phase Locking Value values
% across a range of frequencies for the phase values (x axis) and a range
% of frequencies for the amplitude data (y axis).
% Single electrode analysis.
% Also computes p-values across the PLV surface based on calculated
% surrogate PLV distributions (bootstrapping?)
%
%     SUBID             - Subject's initials (ex- 'ST15'). taken in as a
%                         string
%
%     cond_string       - name of condition to load. taken in as a string.
%                         ex- 'newround_times'
%
%     start_window      - start of window in ms relative to event start
%
%     end_window        - end of window in ms relative to event start
%
%     run_surrogate     - 0 = do not run surrogate analysis.
%                         1 = run surrogate analysis.
%
%     amp_f_array       - 40 frequencies to use for amplitude (optional)
%
%     phase_f_array     - 19 frequencies to use for phase (optional)
%
%
%     Based upon scripts by Ryan Canolty.
%     Modified by Dan Bliss, Sara Szczepanski and Ignacio Saez 11/25/14


%% Initialize
% SUBID = 's01_GP47';
% e_mat = 'gdat_001.mat';

numsurrogate  = 2000; % or more; more is better if you have time to run it
% numsurrogate  = 200; % or more; more is better if you have time to run it

pth_data = ['/Volumes/Data/EcoG_data/gambling_task/',SUBID,'/ordered'];
pth_anal = ['~/Dropbox/Berkeley/Project_EcoG/data/', SUBID, '/'];
pth_img = ['~/Dropbox/Berkeley/Project_EcoG/images/', SUBID, '/'];
load(['/Volumes/Data/EcoG_data/gambling_task/',SUBID,'/signal.mat']); % Takes time but the script is very slow anyways
load(['~/Dropbox/Berkeley/Project_EcoG/data/',SUBID,'/subj_globals.mat']);
e_signal = signal(elec_num,:);
cd(pth_anal);

if ~exist([pth_img 'plv_graphs/'],'dir') % checks if the appropriate subfolder for this study/block has been created in 'analysis'....
    mkdir([pth_img 'plv_graphs/']);
end

if run_surrogate %if we are creating a surrogate distribution
    if ~exist([pth_img 'plv_graphs/' cond_string '/cond_permute_rand_byfreq'],'dir') % checks if the appropriate subfolder for this study/block has been created in 'analysis'....
        mkdir([pth_img 'plv_graphs/' cond_string '/cond_permute_rand_byfreq']);
    end
    
else
    if ~exist([pth_img 'plv_graphs/' cond_string '/cond_permute'],'dir') % checks if the appropriate subfolder for this study/block has been created in 'analysis'....
        mkdir([pth_img 'plv_graphs/' cond_string '/cond_permute']);
    end
    
end


% Read input, and set default values for unsupplied arguments.
for n=1:2:length(varargin)-1
    switch lower(varargin{n})
        case 'amp_f_array'
            amp_f_array = varargin{n+1};
        case 'phase_f_array'
            phase_f_array = varargin{n+1};
    end
end

if ~exist('amp_f_array', 'var')
    amp_f_array = 5:5:200;
end

if ~exist('phase_f_array', 'var')
    phase_f_array = 2:20;
end

% get_signal_parameters, which returns the structure 'sp':
sp = get_signal_parameters('sampling_rate',srate,... % Hz
    'number_points_time_domain',length(e_signal));

n_amp_freqs = length(amp_f_array); %number of frequencies to analyze for amp data
n_phase_freqs = length(phase_f_array);%numer of frequencies to analyze for phase data
plv_matrix_real = nan(n_amp_freqs, n_phase_freqs);%initialize a matrix of plv values
plv_matrix_surr = nan(n_amp_freqs,n_phase_freqs,numsurrogate);%initialize a matrix of plv values- this matrix is as deep as the number of surrogates run.
pvalues_matrix = zeros(n_amp_freqs, n_phase_freqs);%initialize a matrix of p-values for each frequency pair

elec = num2str(elec_num);


%% Analysis
% loops over phase (low) frequencies, then by amp (high) frequencies
% For each combination, it calculates and stores a real PLV value and a
% surrogate PLV value (after randomly circshifting the phase vector)
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
    lf_phase = angle(fs.time_domain);
    clear g.center_frequency
    
    for i_amp = 1:n_amp_freqs %for each frequency of the amplitude data
        
        amp_f = amp_f_array(i_amp); %the particular frequency to filter hf data
        
        g.center_frequency = amp_f;
        g2 = make_chirplet('chirplet_structure', g, 'signal_parameters', sp);
        
        % filter raw signal at high frequency, extract amplitude:
        fs = filter_with_chirplet('raw_signal', e_signal, ...
            'signal_parameters', sp, ...
            'chirplet', g2);
        hf_amp = abs(fs.time_domain);
        % filter high frequency amplitude time-series at low
        % frequency, extract phase:
        fs = filter_with_chirplet('raw_signal', hf_amp, ...
            'signal_parameters', sp, ...
            'chirplet', g1);%filter at low frequency
        hf_phase = angle(fs.time_domain); %extract phase of high frequency amplitude
        
        % this chunks out data for ALL EVENTS in a particular condition AFTER*** filtering
        [chunked_low_phase chunked_high_phase] = get_session_data(SUBID, lf_phase, hf_phase, ...
            cond_string, start_window, end_window);
        
        % Compute cross-frequency phase locking value (PLV).
        %plv = abs(mean(exp(1i*(hf_phase - lf_phase))));
        plv = abs(mean(exp(1i*(chunked_high_phase - chunked_low_phase))));
        
        % Store PLV.
        plv_matrix_real(i_amp, i_phase) = plv;
        
        if run_surrogate %create surrogate distribution of PLV values for each pair of frequencies in the comodulogram
            
            % compute an ensemble of surrogate PLVs to compare to actual value to
            % establish statistical significance:
            iteration = [num2str(i_amp) '_' num2str(i_phase)];
            disp(iteration);
            
            for s = 1:numsurrogate
                
                shift = round(rand*sp.number_points_time_domain); %choose a random timepoint
                surrogate_lf_phase = circshift(chunked_low_phase,[0 shift]); %shift the timecourse of the low frequency phase data at the random cut point
                plv_matrix_surr(i_amp, i_phase, s) = abs(mean(exp(1i*(chunked_high_phase - surrogate_lf_phase))));%surrogate plv value
                clear shift surrogate_lf_phase
                
            end            
            
        end % end if run_surrogate
        
    end % end hf amp for loop
    
end % end lf phase for loop


% if run_surrogate % in case we need to use these for another analysis later
%     save([pth_anal 'plv_graphs/' cond_string '/cond_permute_rand_byfreq/' 'plv_matrix_surr_' elec], 'plv_matrix_surr','plv_matrix_real');
% end

if run_surrogate
    
    % Calculate plv value significance thresholds based upon distribution of
    % plvs for each timepoint; define contour values based upon PLV
    % distribution
    %
    % Returns the PLV values that correspond to certain significance p-values, given a surrogate gamma distribution.
    % These values will then be used as the contour lines on the contourf plot of the real data below.
    %
    %[x_vals, x_vals_all] = pval_from_gamma_dist(surr_matrix)
    %
    [contours, x_vals_all] = pval_from_gamma_dist(plv_matrix_surr);
    
    %save([pth_anal 'plv_graphs/' cond_string '/cond_permute_rand_byfreq/' 'xvals_all_dist_' elec], 'x_vals_all','contours');
end


%% Plot the matrix of PLVs.

% contourf(Z,n) draws a filled contour plot of matrix Z with n contour levels
% contourf(Z,v) draws a filled contour plot of matrix Z with contour lines
%       at the data values specified in the monotonically increasing vector v.
%       The number of contour levels is equal to length(v).


% Note: contourf flips the matrix upside-down, so that the highest
% frequency for amplitude is the top row of the plot (even though
% it's the bottom row of the matrix).

if run_surrogate
    contourf(plv_matrix_real,[0 contours]); %draws plot with specified contour lines
else
    contourf(plv_matrix_real);
end

% Set locations of ticks and their labels.
set(gca, 'XTick', 1:2:19, 'XTickLabel', phase_f_array(1:2:end), ...
    'YTick', 5:5:40, 'YTickLabel', amp_f_array(5:5:end));
title(['Electrode ' elec ' Phase vs. Electrode ' elec ' Amplitude']);
xlabel('frequency for phase (Hz)');
ylabel('frequency for amplitude (Hz)');
h = colorbar;
set(get(h, 'ylabel'), 'string', 'PLV');

if run_surrogate
    print('-dpdf', [pth_img 'plv_graphs/' cond_string '/cond_permute_rand_byfreq/' 'PLV_rand_byfreq_hi_e' elec '_lo_e' elec]);
else
    print('-dpdf', [pth_img 'plv_graphs/' cond_string '/PLV_hi_e' elec '_lo_e' elec]);
    
end

toc

close all

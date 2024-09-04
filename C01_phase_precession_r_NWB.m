
%   Phase precession computation based on spike-time correlation coefficient 

function [] = C01_phase_precession_r_NWB (nwbData, bpfreq)

%% prepare interval data (encoding period)
%=== load data from the interval table
encoding_table = nwbData.intervals.get('encoding_table');

ttls_clip_onsets = encoding_table.start_time.data.load();  % TTL=1
stimCategory = encoding_table.vectordata.get('stimCategory').data.load();  % NB, SB, HB
lfp_data = encoding_table.vectordata.get('lfp').data.load();  % lfp
boundary1_time = encoding_table.vectordata.get('boundary1_time').data.load();  % time of boundary #1
boundary2_time = encoding_table.vectordata.get('boundary2_time').data.load();  % time of boundary #1
boundary3_time = encoding_table.vectordata.get('boundary3_time').data.load();  % time of boundary #1

%=== convert relative boundary times to absolute times by adding start of each video
indxSB = find( stimCategory==1);
indxHB = find( stimCategory==2);

%boundaryTable columns: clipNr/trialNr boundaryType(0,1,2) relativeTime absoluteTime
%== SB
%SB have up to 3 soft boundaries. Add all and prune =0 later.
boundaryTableSB1 = [indxSB ones(length(indxSB),1) boundary1_time(indxSB)];
boundaryTableSB2 = [indxSB ones(length(indxSB),1) boundary2_time(indxSB)];
boundaryTableSB3 = [indxSB ones(length(indxSB),1) boundary3_time(indxSB)];
%SB within HB clips; Each HB clip can have up to two soft boundaries. Add all and prune =0 later
boundaryTableHB_SB1 = [indxHB ones(length(indxHB),1) boundary2_time(indxHB)];
boundaryTableHB_SB2 = [indxHB ones(length(indxHB),1) boundary3_time(indxHB)];
boundaryTableSB_all = [boundaryTableSB1; boundaryTableSB2; boundaryTableSB3; boundaryTableHB_SB1; boundaryTableHB_SB2];
%for SB, prune absoluteTime=0 entires (these dont exist)
indsUse = find( ~isnan(boundaryTableSB_all(:,3)));
%== HB
boundaryTableHB = [indxHB zeros(length(indxHB),1) boundary1_time(indxHB)];
boundaryTable = [boundaryTable; boundaryTableSB_all(indsUse,:) ];
boundaryTable = [boundaryTable; boundaryTableHB];

%% prepare spiking data
baseline_duration = floor(0.5/bin_width_raster);
trial_duration = floor(1/bin_width_raster);

bins_n_raster = baseline_duration + trial_duration;
spks_per_trial = zeros(75, bins_n_raster);
for n_trial = 1:75
    lower_limit = boundary_time_recordered(n_trial)-0.5;
    upper_limit = boundary_time_recordered(n_trial)+ 1;
    spk_indx = (timestampsOfCell > lower_limit) & (timestampsOfCell < upper_limit);
    spk_selected = timestampsOfCell(spk_indx);
    spk_selected_bin = ceil((spk_selected - lower_limit)./(bin_width_raster));
    spks_per_trial(n_trial, spk_selected_bin) = 1;
end

%% Bandpass signal to low frequency bands
cfg = [];
cfg.data            = lfp_data;
cfg.bpfilter        = 'yes';
cfg.bpfreq          = [1 40];
cfg.bpfilttype      = 'fir'; % fFIR filter using MATLAB fir1 function
bp_data_wide = ft_preprocessing(cfg, downsample_data);   
bp_data_wide = bp_data_wide.trial{:};
       
%% Bandpass signal to interested frequency bands
cfg = [];
cfg.bpfilter        = 'yes';
cfg.bpfreq          =  bpfreq;
cfg.bpfilttype      = 'fir'; % fFIR filter using MATLAB fir1 function
bp_data_narrow = ft_preprocessing(cfg, downsample_data);   
bp_data_narrow = bp_data_narrow.trial{:};
    
%% Extract instaneous phase and amplitude from band-passed signals and grouped by conditions
data_phase = angle(hilbert(bp_data_narrow)); % instaneous phase
data_amp = abs(hilbert(bp_data_narrow)); % instaneous amplitude

%% convert time information
bp_data_time = linspace(double(bp_data_narrow.hdr.orig.FirstTimeStamp), double(bp_data_narrow.hdr.orig.LastTimeStamp), bp_data_narrow.sampleinfo(2));
bp_data_time = bp_data_time/10^6;

%% Prep phase and spike information for phase precession calculation
% get TTLs and clip onsets time
clip_onsets = ttls_clip_onsets; % convert to in secs

% extract signals from interested epoches
phase_spike_BAligned = nan(2, downsample_rate*2, length(respMat));
% Dim1: 2 rows, top - theta phases, bottom - binned spike time
% Dim2: number of samples within [-1 1]s window relative to boundaries
% Dim3: number of trials
time_binned = linspace(-1,1, downsample_rate*2+1); % binned time [-1 1]s based on downsample_rate
n_trials = length(boundaryTable);

% extract spiking times, phases relative to [-1, 1]s after boundaries 
for trial_n = 1: n_trials
    boundary_offset = boundaryTable; 
    boundary_time = clip_onsets(trial_n) + boundary_offset;
    window_start_time = boundary_time - 1; 
    window_end_time = boundary_time + 1; 
    try %FIXME later
        phase_spike_BAligned(1,:,trial_n) = data_phase((bp_data_time >= window_start_time) & (bp_data_time <= window_end_time));
    catch
        data_phase_temp = data_phase((bp_data_time >= window_start_time) & (bp_data_time <= window_end_time));
        phase_spike_BAligned(1,:,trial_n) = data_phase_temp(1:downsample_rate*2);
    end
    spike_time_selected = spk_timestamps((spk_timestamps >= window_start_time) & (spk_timestamps <= window_end_time));
    % binned spike time into the same sample size as phase
    edges = linspace(window_start_time, window_end_time, downsample_rate*2+1);
    spike_time_binned = unique(discretize(spike_time_selected, edges)); % spikes happens at the same bin (4ms) are considered as one spike
    phase_spike_BAligned(2,spike_time_binned,trial_n) = time_binned(spike_time_binned); % binned times from [-1 1];       
end 
%% Compute the phase precession using methods from Kempter R. et al, 2012
    time_boundary = reshape(phase_spike_BAligned(2,:,:), [1, n_trials * 2]); 
    phase_boundary = reshape(phase_spike_BAligned(1,:,:), [1, n_trials * 2]);         
    [cc, p, slope, offset, MRL] = cl_corr(time_boundary, phase_boundary);
    
 %% plot phase precession results
    figure('rend','painters','pos',[10 10 1300 600]);
    scatter(time_boundary, phase_boundary, 20, [0 0 0],'filled'); hold on
    ylim([-pi pi]);
    xlim([0 3*pi])
    ylabel('Theta phase (o)')
    set(gca, 'XTick', [0 3*pi], 'XTickLabel',[0 3*pi],  'YTick',[-pi 0 pi], 'YTickLabel', rad2deg([-pi 0 pi]), ....
             'FontSize', 20, 'FontWeight', 'bold', 'box', 'on', 'LineWidth', 2)
    title(['r = ', num2str(cc), ',p = ', num2str(p), ', slope= ', num2str(slope), ', offset= ', num2str(offset), ', MRL= ',num2str(MRL)]);
end
    
















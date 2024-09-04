
%   Theta bouts detection
%   Modified by JZ based on "Bicycle" analyses suite: https://github.com/bycycle-tools/bycycle
%
% Parameters settings
% 
%   cycleLength= 3;      % least cycle numbers
%   amplitude = 0.3;           % amplitude threshold, 0.3
%   period = 0.5;        % period threshold, 0.5
%   mono = 0.6;          % Monocity threshold 0.6
%   bpfreq = [3 4];
%   similar parameters as used for motor beta recording in Cole et al., 2019

function [] = C01_thetaBouts_detection_NWB (cyclelength, amplitude, period, mono, bpfreq)

% Parameters settings
params.CycleLength= cyclelength;      % least cycle numbers
params.Amp = amplitude;           % amplitude threshold, 0.3
params.Period = period;        % period threshold, 0.5
params.Mono = mono;          % Monocity threshold 0.6
params.bpfreq = bpfreq;

%% prepare interval data (encoding period)
% === load data from the interval table
encoding_table = nwbData.intervals.get('encoding_table');
ttls_clip_onsets = encoding_table.start_time.data.load();  % TTL=1
lfp_data = encoding_table.vectordata.get('lfp').data.load();  

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
   
%% extract the zero crossings
ZeroCrosser = ZeroX(1:length(bp_data_narrow), bp_data_narrow);
ZeroCrosser = floor(ZeroCrosser); % convert to integers (indx)
ZeroCrosser = unique(ZeroCrosser);% remove identical adjacent zero crossing points;
    
%% find troughs and peaks
peaks = [];
troughs = [];
for zero_n = 1:length(ZeroCrosser)-1
    data_win = bp_data_wide(ZeroCrosser(zero_n):ZeroCrosser(zero_n+1));
    if sum(data_win) > 0
       [~, peak_indx] = max(data_win);
       peak_indx = peak_indx + ZeroCrosser(zero_n) -1;
       peaks = [peaks, peak_indx];
    elseif sum(data_win) < 0
       [~, trough_indx] = min(data_win);
       trough_indx = trough_indx + ZeroCrosser(zero_n) -1;
       troughs = [troughs, trough_indx];
    end        
end 
    
% for simplicity, ensure each cycle starts with trough, zero, peak and zero, trough 
% create arrays for saying critical points for good cycles
good_cycle_n = 0;
% Dim1: peak_time vs trough_time vs zeroCrosser_time
% Dim2: first vs second
    
for trough_n = 2:length(troughs)-1 % skip the first cycle, otherwise the Peak-trough symmetry cannot be computed
    % Ensure every cycle includes one peak and two zeros,
    trough_first = troughs(trough_n);
    trough_second = troughs(trough_n + 1);
    [~, peak_pick] = find(peaks > trough_first & peaks < trough_second);
    [~, zeros_pick] = find(ZeroCrosser > trough_first & ZeroCrosser < trough_second);
    % compute oscillatory features:
    if length(peak_pick) == 1 && length(zeros_pick) == 2
       good_cycle_n = good_cycle_n +1;
       % Amplitude: the avg value from the two adjacent troughs around one peak
       Amp1 = bp_data_wide(peaks(peak_pick)) - bp_data_wide(trough_first);
       Amp2 = bp_data_wide(peaks(peak_pick)) - bp_data_wide(trough_second);
       AmpArray(trough_n)= mean([Amp1,Amp2]);
       % From the paper: the amplitude consistency of a cycle is quanti?ed as the relative difference in 
       % the rise and decay voltage (e.g., 0.5 corresponds to the change in voltage in one
       % ?ank being 2 times bigger than the other, 1.0 corresponds to the rise and decay 
       % ?anks having equal changes in voltage, etc.). The minimum value is taken after this
       % measure is computed for each pair of adjacent rise and decay ?anks that includes one of the current cycle?s ?anks.
       
       % I think such definition for amplitude consistency is
       % problematic, as the ratio can be flipping across cycles,
       % for example, one cycle rise = 2*decay but next cycle decay = 2*rise;
       if Amp1/Amp2 <= 1
          AmpConsistency(trough_n) = Amp1/Amp2;
       else
          AmpConsistency(trough_n) = Amp2/Amp1; 
       end
       % Period: the time between two consecutive troughs, in data points, check fs to convert to seconds
       period_rise = peaks(peak_pick) - trough_first;
       period_decay = trough_second - peaks(peak_pick);
       PeriodArray(trough_n)= period_rise + period_decay;
       PeriodConsistency(trough_n) = period_rise/period_decay;
       % Rise-decay symmetry: ratio of period_rise in each cycle
       rdsymArray(trough_n) = period_rise/(period_rise + period_decay);
       % Peak-trough symmetry: ratio of period_trough from perious cycle
       period_trough = troughs(trough_n) - troughs(trough_n - 1);
       period_peak = troughs(trough_n + 1) - troughs(trough_n);
       ptsymArray(trough_n) = period_trough/(period_trough + period_peak);
       % Monotonicity: the fraction of instantaneous voltage changes (difference between consecutive samples)
       % that are positive during the rise phase and negative during the decay phase
       % a sine wave would have a perfect 1
       % 0.8 corresponds to 20 of the voltage time series going in the opposite direction of the current flank
       mono_rise = sum(diff(bp_data_wide(trough_first:peaks(peak_pick))) > 0)/(peaks(peak_pick) - trough_first);
       mono_decay = sum(diff(bp_data_wide(peaks(peak_pick):trough_second)) < 0)/(trough_second - peaks(peak_pick));
       MonoConsistency(trough_n) = mono_rise/mono_decay;
       % cycle times
       %PeriodTimes(trough_n,1)= bp_data_wide(troughs(trough_n)); 
       %PeriodTimes(trough_n,2)= bp_data_wide(troughs(trough_n+1));
       good_cycles_points(1, 1, good_cycle_n) = peaks(peak_pick);
       good_cycles_points(2, 1, good_cycle_n) = trough_first;
       good_cycles_points(2, 2, good_cycle_n) = trough_second;
       good_cycles_points(3, 1, good_cycle_n) = ZeroCrosser(zeros_pick(1));
       good_cycles_points(3, 2, good_cycle_n) = ZeroCrosser(zeros_pick(2));
    end
end    

% interpolating phase information,
% trough_first - pi, zeros_first = -pi/2, peaks - 0 degree; zeros_second pi/2; trough_second: pi
PhaseArray=nan(1, length(bp_data_wide));
for cycle_n=1:good_cycle_n-1
    Phases_rise=linspace(-pi,0, (good_cycles_points(1,1,cycle_n)- good_cycles_points(2,1,cycle_n)+1));    
    Phases_decay=linspace(0,pi, (good_cycles_points(2,2,cycle_n)- good_cycles_points(1,1,cycle_n)+1));    
    PhaseArray(good_cycles_points(2,1,cycle_n):good_cycles_points(1,1,cycle_n)-1)= Phases_rise(1:end-1);
    PhaseArray(good_cycles_points(1,1,cycle_n):good_cycles_points(2,2,cycle_n))= Phases_decay;
end
% detect cycles with consistent amplitude, period, and monotonicity above the threshold
amp_pass_indx = find(AmpConsistency >= params.Amp);
period_pass_indx = find(PeriodConsistency >= params.Period);
mono_pass_indx = find(MonoConsistency >= params.Mono);
semi_pass_indx = intersect(intersect(amp_pass_indx, period_pass_indx), mono_pass_indx);
% detect oscillations longer than the minimum cycle length
semi_pass_indx_diff = diff(semi_pass_indx);
semi_pass_indx_diff(semi_pass_indx_diff ~= 1) = 0;
c = movprod(semi_pass_indx_diff, 2);
pass_indx= find(c == 1);
final_pass_indx = unique([pass_indx-1, pass_indx, pass_indx+1]);
final_pass_indx(final_pass_indx == 0) = [];
osc_trough_selected = squeeze(good_cycles_points(2, 1, final_pass_indx));
    
%% plotting detected theta bouts
   for trough_n = 1:10:length(final_pass_indx) % FIXME later
       start_sample = good_cycles_points(2, 1, final_pass_indx(trough_n))- downsample_rate*1;
       end_sample = good_cycles_points(2, 1, final_pass_indx(trough_n)) + downsample_rate*3; % plot 5s
       n_cycles = sum(osc_trough_selected > start_sample & osc_trough_selected < end_sample);
       fig = figure ('renderer', 'painters', 'position', [10 10 2000 800], 'visible', 'on');
       subplot(3,1,1)
       plot(downsample_data.trial{1,1}(start_sample:end_sample), 'Color', [0.5 0.5 0.5], 'LineWidth', 3); hold on % raw signal, 
       plot(bp_data_wide(start_sample:end_sample),'Color', [0 0 0], 'LineWidth', 1.5); % band-pass filtered signals
       for cycle_n = 1:n_cycles
           time = good_cycles_points(2, 1, final_pass_indx(trough_n)+cycle_n-1): good_cycles_points(2, 2, final_pass_indx(trough_n)+cycle_n-1); 
           plot(time-start_sample+1, bp_data_wide(time), 'r', 'LineWidth', 2); % detected oscillations
       end
       xlim([0 4]*downsample_rate);
       ylabel('Amplitude (uV)');
       legend({'raw signal', 'Low-passed', 'Detected oscillations'});
       title(['Cycle#', num2str(final_pass_indx(trough_n)), '-', num2str(final_pass_indx(trough_n+n_cycles-1)), ...
               '  Amp=', num2str(params.Amp), '  Period=', num2str(params.Period),...
               '  Mono=', num2str(params.Mono), '  CycleLength>=', num2str(params.CycleLength) ...
               '  bp=', num2str(bpfreq)]);
       set(gca, 'xtick', [0:4].*downsample_rate, 'xtickLabel', 0:4, 'box', 'on', 'LineWidth', 1.5 , 'fontSize', 20);
       subplot(3,1,2)
       plot(bp_data_narrow(start_sample:end_sample),'Color', [0 0 0], 'LineWidth', 1.5); % band-pass filtered signals
       set(gca, 'xtick', [0:4].*downsample_rate, 'xtickLabel', 0:4, 'box', 'on', 'LineWidth', 1.5 , 'fontSize', 20);
       xlim([0 4]*downsample_rate);
       ylabel('Amplitude (uV)'); 
       legend('Theta-filtered')
       subplot(3,1,3)
       plot(data_phase(start_sample:end_sample),'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); hold on % Hilbert-way
       for cycle_n = 1:n_cycles
           time = good_cycles_points(2, 1, final_pass_indx(trough_n)+cycle_n-1): good_cycles_points(2, 2, final_pass_indx(trough_n)+cycle_n-1); 
           plot(time-start_sample+1, PhaseArray(time), 'r', 'LineWidth', 1.5); % detected oscillations
       end
       set(gca, 'xtick', [0:4].*downsample_rate, 'xtickLabel', 0:4, 'box', 'on', 'LineWidth', 1.5 , 'fontSize', 20);
       xlim([0 4]*downsample_rate);
       xlabel('Time');
       legend({'Hilbert-phase', 'Cycle-phase'})
       chan_indx = sscanf(chan_list(chan_n).name, 'CSC%d.nsc');
       figure_name = [figure_dir, '/', subject_id, '_A_ss', num2str(chan_indx) , '_cycle', num2str(final_pass_indx(trough_n)), 'to', num2str(final_pass_indx(trough_n+n_cycles-1)), '.png'];
       saveas(fig, figure_name); % FIXME later
       close(fig);
   end  
end
    
















%
%
% B01_raster_psth_timeDiscrim_SAligned_BSep_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag)
%
%History
%12/7/21: This file is forked from original B01_raster_psth_timeDiscrim_BAligned_BSep.m by JZ.
%

function [] = B01_raster_psth_timeDiscrim_RAligned_BSep_NWB(nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag)

%% prepare interval data (time discrimination period)
%=== load data from the interval table
timediscrimination_table = nwbData.intervals.get('timediscrimination_table');

ttls_img_offsets = timediscrimination_table.stop_time.data.load();  % TTL=2
response = timediscrimination_table.vectordata.get('RT').data.load() + ttls_img_offsets;
stimCategory = timediscrimination_table.vectordata.get('boundary_type').data.load();  % NB, SB, HB
accuracy = timediscrimination_table.vectordata.get('accuracy').data.load();  % Correct, incorrect


%=== convert relative boundary times to absolute times by adding start of each video
indxNB_CR = ( stimCategory==1 & accuracy==2);
indxSB_CR = ( stimCategory==2 & accuracy==2);
indxHB_CR = ( stimCategory==3 & accuracy==2);
indxNB_IC = ( stimCategory==1 & accuracy==1);
indxSB_IC = ( stimCategory==2 & accuracy==1);
indxHB_IC = ( stimCategory==3 & accuracy==1);

%% prepare spiking data
baseline_duration = floor(1/bin_width_raster);
trial_duration = floor(0.5/bin_width_raster);

bins_n_raster = baseline_duration + trial_duration;
spks_per_trial = zeros(180, bins_n_raster);
for n_trial = 1:180
    lower_limit = response(n_trial)-1;
    upper_limit = response(n_trial)+0.5;
    spk_indx = (timestampsOfCell > lower_limit) & (timestampsOfCell < upper_limit);
    spk_selected = timestampsOfCell(spk_indx);
    spk_selected_bin = ceil((spk_selected - lower_limit)./(bin_width_raster));
    spks_per_trial(n_trial, spk_selected_bin) = 1;
end

%% plot

subjectID = nwbData.general_session_id; %CBID

cellLabelStr = ['NWB ' subjectID '-' num2str(channelid) '-' num2str(cellNr) '-' num2str(brainAreaOfCell)];

[xpoints_NB_CR,ypoints_NB_CR] = find(spks_per_trial(indxNB_CR,:) == 1); % NB
[xpoints_SB_CR,ypoints_SB_CR] = find(spks_per_trial(indxSB_CR,:) == 1); % SB
[xpoints_HB_CR,ypoints_HB_CR] = find(spks_per_trial(indxHB_CR,:) == 1); % HB
[xpoints_NB_IC,ypoints_NB_IC] = find(spks_per_trial(indxNB_IC,:) == 1); % NB
[xpoints_SB_IC,ypoints_SB_IC] = find(spks_per_trial(indxSB_IC,:) == 1); % SB
[xpoints_HB_IC,ypoints_HB_IC] = find(spks_per_trial(indxHB_IC,:) == 1); % HB

HB_color = [255,99,71]./255;
SB_color = [100,149,237]./255;
NB_color = [102 204 0]./255;
fig = figure('rend','painters','pos',[10 10 400 700], 'visible', plot_flag);
h_raster = subplot(2,1,1);
ax_raster = get(h_raster,'Position');
ax_raster(4) = ax_raster(4)+0.1;
ax_raster(2) = ax_raster(2)-0.12;
set(h_raster, 'Position', ax_raster);
scatter(ypoints_NB_CR, xpoints_NB_CR-min(xpoints_NB_CR)+121, 4, 'MarkerEdgeColor', NB_color, 'MarkerFaceColor', NB_color, 'LineWidth', 1); hold on
scatter(ypoints_NB_IC, xpoints_NB_IC-min(xpoints_NB_IC)+sum(indxNB_CR)+121, 8, 'MarkerEdgeColor', NB_color, 'MarkerFaceColor', [1 1 1], 'LineWidth', 1);
scatter(ypoints_SB_CR, xpoints_SB_CR-min(xpoints_SB_CR)+61, 4, 'MarkerEdgeColor', SB_color, 'MarkerFaceColor', SB_color, 'LineWidth', 1);
scatter(ypoints_SB_IC, xpoints_SB_IC-min(xpoints_SB_IC)+sum(indxSB_CR)+61, 8, 'MarkerEdgeColor', SB_color, 'MarkerFaceColor', [1 1 1], 'LineWidth', 1);
scatter(ypoints_HB_CR, xpoints_HB_CR-min(xpoints_HB_CR)+1, 4, 'MarkerEdgeColor', HB_color, 'MarkerFaceColor', HB_color, 'LineWidth', 1);
scatter(ypoints_HB_IC, xpoints_HB_IC-min(xpoints_HB_IC)+sum(indxHB_CR)+1, 8, 'MarkerEdgeColor', HB_color, 'MarkerFaceColor', [1 1 1], 'LineWidth', 1);
plot([1/bin_width_raster 1/bin_width_raster], [0 180],'k--','LineWidth',2); % the trial onset
xlim([0 baseline_duration + trial_duration])
ylim([0 180])
xlabel('Time (seconds)')
ylabel('Trial number')
set(gca, 'LineWidth', 1.5, 'XTick', 11/bin_width_raster*[0 1 1.5], 'XTickLabel', [-1 0 0.5],...
   'fontsize', 15, 'fontweight', 'bold','box','on')    
% generate firing rate plots, 200ms window
n_data = (bin_width_fr/bin_width_raster); % number of data per window
fr_per_bin_HB_CR = smoothdata(spks_per_trial(indxHB_CR,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_SB_CR = smoothdata(spks_per_trial(indxSB_CR,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_NB_CR = smoothdata(spks_per_trial(indxNB_CR,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_HB_IC = smoothdata(spks_per_trial(indxHB_IC,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_SB_IC = smoothdata(spks_per_trial(indxSB_IC,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_NB_IC = smoothdata(spks_per_trial(indxNB_IC,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_avg_HB_CR = mean(fr_per_bin_HB_CR,1);
fr_per_bin_avg_SB_CR = mean(fr_per_bin_SB_CR,1);
fr_per_bin_avg_NB_CR = mean(fr_per_bin_NB_CR,1);
fr_per_bin_avg_HB_IC = mean(fr_per_bin_HB_IC,1);
fr_per_bin_avg_SB_IC = mean(fr_per_bin_SB_IC,1);
fr_per_bin_avg_NB_IC = mean(fr_per_bin_NB_IC,1);
h_psth = subplot(2,1,2);
ax_psth = get(h_psth,'Position');
ax_psth(4) = ax_psth(4)-0.1;
ax_psth(2) = ax_psth(2)-0.01;
set(h_psth, 'Position', ax_psth);
boundedline(1:bins_n_raster,fr_per_bin_avg_HB_CR , std(fr_per_bin_HB_CR,1)./sqrt(sum(indxHB_CR)),'cmap',HB_color,'alpha'); hold on
boundedline(1:bins_n_raster,fr_per_bin_avg_HB_IC , std(fr_per_bin_HB_IC,1)./sqrt(sum(indxHB_IC)), '-.', 'cmap',HB_color,'alpha');
boundedline(1:bins_n_raster,fr_per_bin_avg_SB_CR , std(fr_per_bin_SB_CR,1)./sqrt(sum(indxSB_CR)),'cmap',SB_color,'alpha');
boundedline(1:bins_n_raster,fr_per_bin_avg_SB_IC , std(fr_per_bin_HB_IC,1)./sqrt(sum(indxSB_IC)), '-.', 'cmap',SB_color,'alpha');
boundedline(1:bins_n_raster,fr_per_bin_avg_NB_CR , std(fr_per_bin_NB_CR,1)./sqrt(sum(indxNB_CR)),'cmap',NB_color,'alpha');
boundedline(1:bins_n_raster,fr_per_bin_avg_NB_IC , std(fr_per_bin_NB_IC,1)./sqrt(sum(indxNB_IC)), '-.', 'cmap',NB_color,'alpha');
plot([1/bin_width_raster 1/bin_width_raster], [0 180],'k--','LineWidth',2); % the trial onset
xlim([0 baseline_duration + trial_duration])
ylim([0 max([fr_per_bin_avg_HB_CR, fr_per_bin_avg_SB_CR, fr_per_bin_avg_NB_CR fr_per_bin_avg_HB_IC, fr_per_bin_avg_SB_IC, fr_per_bin_avg_NB_IC]+2)])
xlabel('Time (seconds)')
ylabel('Firing rate (Hz)')
set(gca, 'LineWidth', 1.5, 'XTick',  1/bin_width_raster*[0 1 1.5], 'XTickLabel', [-1 0 0.5],...
        'fontsize', 15, 'fontweight', 'bold','box','on')

title(cellLabelStr, 'Interpreter', 'none');
figure_dir = [pwd, '/NWB/exported/Figures/CBID', num2str(subjectID), '/RAligned_timeDiscrim/'];
if ~exist(figure_dir)
    mkdir(figure_dir)
end

figure_name = ['Raster_psth_RAligned', cellLabelStr, '_timeDiscrim.png'];
print(fig,'-dpng',[figure_dir,'/', figure_name]);
close(fig)


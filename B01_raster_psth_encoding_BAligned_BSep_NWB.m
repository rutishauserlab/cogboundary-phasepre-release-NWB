%
%
% B01_raster_psth_encoding_BAligned_BSep_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag)
%
%History
%12/2/21: This file is forked from original B01_raster_psth_encoding_BAligned_BSep.m by JZ.
%

function [] = B01_raster_psth_encoding_BAligned_BSep_NWB(nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag)

%% prepare interval data (encoding period)
%=== load data from the interval table
encoding_table = nwbData.intervals.get('encoding_table');

ttls_clip_onsets = encoding_table.start_time.data.load();  % TTL=1
stimCategory = encoding_table.vectordata.get('stimCategory').data.load();  % NB, SB, HB
boundary1_time = encoding_table.vectordata.get('boundary1_time').data.load();  % time of boundary #1
boundary2_time = encoding_table.vectordata.get('boundary2_time').data.load();  % time of boundary #1
boundary3_time = encoding_table.vectordata.get('boundary3_time').data.load();  % time of boundary #1

%=== convert relative boundary times to absolute times by adding start of each video
indxNB = find( stimCategory==0);
indxSB = find( stimCategory==1);
indxHB = find( stimCategory==2);

%boundaryTable columns: clipNr/trialNr boundaryType(0,1,2) relativeTime absoluteTime
%== NB
boundaryTableNB = [indxNB zeros(length(indxNB),1) boundary1_time(indxNB)];

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


boundaryTable = [boundaryTableNB];
boundaryTable = [boundaryTable; boundaryTableSB_all(indsUse,:) ];
boundaryTable = [boundaryTable; boundaryTableHB];

% Modified by JZ, boundary time and spike time has already aligned to the
% task onset, no need to correct for absolute time
% %=== convert from relative to absolute time 
% for k=1:size(boundaryTable,1)
%     boundaryTable(k,4) = ttls_clip_onsets(boundaryTable(k,1)) + boundaryTable(k,3);
% end

%boundary_time_recordered  contains trials re-ordered - NB, SB, HB
% boundary_time_recordered = boundaryTable(:,4)';

boundary_time_recordered = boundaryTable(:,3);

%% prepare spiking data
baseline_duration = floor(0.5/bin_width_raster);
trial_duration = floor(1/bin_width_raster);

bins_n_raster = baseline_duration + trial_duration;
spks_per_trial = zeros(135, bins_n_raster);
for n_trial = 1:135
    lower_limit = boundary_time_recordered(n_trial)-0.5;
    upper_limit = boundary_time_recordered(n_trial)+ 1;
    spk_indx = (timestampsOfCell > lower_limit) & (timestampsOfCell < upper_limit);
    spk_selected = timestampsOfCell(spk_indx);
    spk_selected_bin = ceil((spk_selected - lower_limit)./(bin_width_raster));
    spks_per_trial(n_trial, spk_selected_bin) = 1;
end

%% plot

subjectID = nwbData.general_session_id; %CBID

cellLabelStr = ['NWB ' subjectID '-' num2str(channelid) '-' num2str(cellNr) '-' num2str(brainAreaOfCell)];


[xpoints_NB,ypoints_NB] = find(spks_per_trial(1:30,:) == 1); % NB
[xpoints_SB,ypoints_SB] = find(spks_per_trial(31:105,:) == 1); % SB
[xpoints_HB,ypoints_HB] = find(spks_per_trial(106:135,:) == 1); % HB

HB_color = [255,99,71]./255;
SB_color = [100,149,237]./255;
NB_color = [102 204 0]./255;
fig = figure('rend','painters','pos',[10 10 400 700], 'visible', plot_flag);
h_raster = subplot(2,1,1);
ax_raster = get(h_raster,'Position');
ax_raster(4) = ax_raster(4)+0.1;
ax_raster(2) = ax_raster(2)-0.12;
set(h_raster, 'Position', ax_raster);
scatter(ypoints_NB, xpoints_NB-min(xpoints_NB)+105, 1, 'MarkerEdgeColor', NB_color, 'MarkerFaceColor', NB_color, 'LineWidth', 1); hold on
scatter(ypoints_SB, xpoints_SB-min(xpoints_SB)+31, 1, 'MarkerEdgeColor', SB_color, 'MarkerFaceColor', SB_color, 'LineWidth', 1);
scatter(ypoints_HB, xpoints_HB-min(xpoints_HB)+1, 1, 'MarkerEdgeColor', HB_color, 'MarkerFaceColor', HB_color, 'LineWidth', 1);
plot([0.5/bin_width_raster 0.5/bin_width_raster], [0 135], 'k--', 'LineWidth', 2);
xlim([0 baseline_duration + trial_duration])
ylim([0 135])
xlabel('Time (seconds)')
ylabel('Trial number')
set(gca, 'LineWidth', 1.5, 'XTick', 1/bin_width_raster*[0 0.5 1.5], 'XTickLabel', [-0.5 0 1],...
    'fontsize', 15, 'fontweight', 'bold','box','on')
% generate firing rate plots, 200ms window
n_data = (bin_width_fr/bin_width_raster); % number of data per window
fr_per_bin_HB = smoothdata(spks_per_trial(106:135,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_SB = smoothdata(spks_per_trial(31:105,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_NB = smoothdata(spks_per_trial(1:30,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_avg_HB = mean(fr_per_bin_HB,1);
fr_per_bin_avg_SB = mean(fr_per_bin_SB,1);
fr_per_bin_avg_NB = mean(fr_per_bin_NB,1);
h_psth = subplot(2,1,2);
ax_psth = get(h_psth,'Position');
ax_psth(4) = ax_psth(4)-0.1;
ax_psth(2) = ax_psth(2)-0.01;
set(h_psth, 'Position', ax_psth);
boundedline(1:bins_n_raster,fr_per_bin_avg_HB , std(fr_per_bin_HB,1)./sqrt(30),'cmap',HB_color,'alpha'); hold on
boundedline(1:bins_n_raster,fr_per_bin_avg_SB , std(fr_per_bin_SB,1)./sqrt(75),'cmap',SB_color,'alpha');
boundedline(1:bins_n_raster,fr_per_bin_avg_NB , std(fr_per_bin_NB,1)./sqrt(30),'cmap',NB_color,'alpha');
plot([0.5/bin_width_raster 0.5/bin_width_raster], [0 135],'k--','LineWidth',2); % the trial onset
xlim([0 baseline_duration + trial_duration])
ylim([0 max([fr_per_bin_avg_HB, fr_per_bin_avg_SB, fr_per_bin_avg_NB]+2)])
xlabel('Time (seconds)')
ylabel('Firing rate (Hz)')
set(gca, 'LineWidth', 1.5, 'XTick', 1/bin_width_raster*[0 0.5 1.5], 'XTickLabel', [-0.5 0 1],...
    'fontsize', 15, 'fontweight', 'bold','box','on')

title(cellLabelStr, 'Interpreter', 'none');
figure_dir = [pwd, '/NWB/exported/Figures/CBID', num2str(subjectID), '/BAligned_encoding/'];
if ~exist(figure_dir)
    mkdir(figure_dir)
end

figure_name = ['Raster_psth_BAligned', cellLabelStr, '_encoding.png'];
print(fig,'-dpng',[figure_dir,'/', figure_name]);
close(fig)


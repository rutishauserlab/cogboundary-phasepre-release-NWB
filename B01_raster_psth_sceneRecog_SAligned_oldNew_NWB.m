%
%
% B01_raster_psth_sceneRecog_BAligned_BSep_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag)
%
%History
%12/2/21: This file is forked from original B01_raster_psth_sceneRecog_BAligned_BSep.m by JZ.
%

function [] = B01_raster_psth_sceneRecog_SAligned_oldNew_NWB(nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag)

%% prepare interval data (scene recognition period)
%=== load data from the interval table
recognition_table = nwbData.intervals.get('recognition_table');

ttls_img_onsets = recognition_table.start_time.data.load();  % TTL=1

accuracy = recognition_table.vectordata.get('accuracy').data.load();  % Correct, incorrect
old_new  = recognition_table.vectordata.get('old_new').data.load();  % Correct, incorrect

plotMode = 2; % 1 correct/incorrect,  2 old/new
if plotMode==1
    % == accuracy. plot correct vs incorrect response
    indxGrp1 = ( accuracy==2);
    indxGrp2 = ( accuracy==1);
    grpLabels = {'Corr','Incorr'};
else
    % == old vs new  (correct only)
    indxGrp1 = ( old_new==2 & accuracy==2);
    indxGrp2 = ( old_new==1 & accuracy==2);
    grpLabels = {'Old','New'};
end

%% prepare spiking data
baseline_duration = floor(0.5/bin_width_raster);
trial_duration = floor(3.5/bin_width_raster);

bins_n_raster = baseline_duration + trial_duration;
spks_per_trial = zeros(180, bins_n_raster);
for n_trial = 1:180
    lower_limit = ttls_img_onsets(n_trial)-0.5;
    upper_limit = ttls_img_onsets(n_trial)+ 3.5;
    spk_indx = (timestampsOfCell > lower_limit) & (timestampsOfCell < upper_limit);
    spk_selected = timestampsOfCell(spk_indx);
    spk_selected_bin = ceil((spk_selected - lower_limit)./(bin_width_raster));
    spks_per_trial(n_trial, spk_selected_bin) = 1;
end

%% plot

subjectID = nwbData.general_session_id; %CBID

cellLabelStr = ['NWB ' subjectID '-' num2str(channelid) '-' num2str(cellNr) '-' num2str(brainAreaOfCell)];

[xpoints_CR,ypoints_CR] = find(spks_per_trial(indxGrp1,:) == 1); % NB
[xpoints_IC,ypoints_IC] = find(spks_per_trial(indxGrp2,:) == 1); % SB

fig = figure('rend','painters','pos',[10 10 1200 700], 'visible', plot_flag);
h_raster = subplot(2,1,1);
ax_raster = get(h_raster,'Position');
ax_raster(4) = ax_raster(4)+0.1;
ax_raster(2) = ax_raster(2)-0.12;
set(h_raster, 'Position', ax_raster);
scatter(ypoints_CR, xpoints_CR-min(xpoints_CR)+1, 4, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1); hold on
scatter(ypoints_IC, xpoints_IC-min(xpoints_IC)+1+sum(indxGrp1), 8, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 1], 'LineWidth', 1);
plot([0.5/bin_width_raster 0.5/bin_width_raster], [0 180],'k--','LineWidth',2); % the image onset
plot([2/bin_width_raster 2/bin_width_raster], [0 180],'k--','LineWidth',2); % the image offset
xlim([0 baseline_duration + trial_duration])
ylim([0 180])
xlabel('Time (seconds)')
ylabel('Trial number')
set(gca, 'LineWidth', 1.5, 'XTick', 1/bin_width_raster*[0 0.5 2 4], 'XTickLabel', [-0.5 0 1.5 3.5],...
   'fontsize', 15, 'fontweight', 'bold','box','on')    
% generate firing rate plots, 200ms window
n_data = (bin_width_fr/bin_width_raster); % number of data per window
fr_per_bin_IC = smoothdata(spks_per_trial(indxGrp2,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_CR = smoothdata(spks_per_trial(indxGrp1,:),2, 'gaussian',n_data)*n_data;  % convert it to Hz
fr_per_bin_avg_IC = mean(fr_per_bin_IC,1);
fr_per_bin_avg_CR = mean(fr_per_bin_CR,1);
h_psth = subplot(2,1,2);
ax_psth = get(h_psth,'Position');
ax_psth(4) = ax_psth(4)-0.1;
ax_psth(2) = ax_psth(2)-0.01;
set(h_psth, 'Position', ax_psth);

hl_Grp1 = boundedline(1:bins_n_raster,fr_per_bin_avg_IC , std(fr_per_bin_IC,1)./sqrt(sum(indxGrp2)),'-.', 'cmap',[0 0 0],'alpha'); hold on
hl_Grp2 = boundedline(1:bins_n_raster,fr_per_bin_avg_CR , std(fr_per_bin_CR,1)./sqrt(sum(indxGrp1)),'cmap',[0 0 0],'alpha');

plot([0.5/bin_width_raster 0.5/bin_width_raster], [0 180],'k--','LineWidth',2); % the trial onset
plot([2/bin_width_raster 2/bin_width_raster], [0 180],'k--','LineWidth',2); % the trial offset
xlim([0 baseline_duration + trial_duration])
ylim([0 max([fr_per_bin_avg_IC, fr_per_bin_avg_CR]+2)])
xlabel('Time (seconds)')
ylabel('Firing rate (Hz)')
set(gca, 'LineWidth', 1.5, 'XTick', 1/bin_width_raster*[0 0.5 2 4], 'XTickLabel', [-0.5 0 1.5 3.5],...
        'fontsize', 15, 'fontweight', 'bold','box','on')

title(cellLabelStr, 'Interpreter', 'none');
figure_dir = [pwd, '/NWB/exported/Figures/CBID', num2str(subjectID), '/SAligned_sceneRecog_oldNew/'];
if ~exist(figure_dir)
    mkdir(figure_dir)
end

legend([hl_Grp1(1) hl_Grp2(1)], grpLabels);

figure_name = ['Raster_psth_SAligned', cellLabelStr, '_sceneRecog_oldNew.png'];
print(fig,'-dpng',[figure_dir,'/', figure_name]);
close(fig)


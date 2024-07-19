%
% Read and plot cognitive boundary data based on NWB-formatted files
%
% Pre-requisite: NWB routines are installed and in path

%% Section 1 -- Set paths and initialize NWB API

basepathData = [pwd, '/exported_upload/'];  % NWB data files

basePathNWBCode = [pwd, '/matnwb-2.4.0.0/'];
basePathCode = [pwd, '/NWB_export_JZ_new/'];

%== Initialize NWB
if exist([basePathNWBCode, filesep, '+types', filesep, '+core', filesep, 'NWBFile.m'])
     disp(['generateCore() already initialized...']) %only need to do once
else 
    cd([basePathNWBCode])
    generateCore();
    disp(['generateCore() initialized...'])
end 

addpath(basePathNWBCode); 
addpath([basePathCode filesep 'export'  ]); 
addpath([basePathCode filesep 'export' filesep 'helpers' ]); 
addpath([basePathCode filesep 'analysisNWB' ]);

%% open NWB file for one session and import spiking data needed for plotting
%which session to analyze
for n = 1:19 % changed by JZ on Dec 1, 2021
% for n = 6  
% fName_in = 'CBID8.nwb';
fName_in = ['CBID', num2str(n), '.nwb']; % changed by JZ on Dec 1, 2021

% import file
nwbData = nwbRead([basepathData filesep fName_in]);

% load neural data
all_spike_data = nwbData.units.spike_times.data.load();
spike_data_indexes = nwbData.units.spike_times_index.data.load();
channel_ids_index = nwbData.general_extracellular_ephys_electrodes.vectordata.get('origChannel').data.load();
cell_electrodes = nwbData.units.electrodes.data.load();
brain_areas_index = nwbData.general_extracellular_ephys_electrodes.vectordata.get('location').data.load;
clusterIDs = nwbData.units.id.data.load;   %all available id's (clusters)
elec_index = [nwbData.units.electrodes.data.load()]+1;
channel_ids = channel_ids_index(elec_index); 
brain_areas = brain_areas_index(elec_index,:); 

bin_width_raster = 0.002; % 10ms per bin
bin_width_fr = 0.2; % 200ms per bin
plot_flag = 'on';

    %loop over all cells in this file
    for i = 1:length(clusterIDs)
        %== locate spikes that belong to this cluster
        timestampsOfCell = nwb_read_unit(nwbData.units.spike_times_index, nwbData.units.spike_times, clusterIDs(i)+1 );

        channelid =  channel_ids(i);
        % cellNr = cell_electrodes(i);
        cellNr = clusterIDs(i); % made by JZ on Dec 1, 2021
        brainAreaOfCell = brain_areas(i, :);
        currentCell = cellNr; 

        %== plot this cell
        B01_raster_psth_encoding_BAligned_BSep_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag);
        
        %B01_raster_psth_encoding_SAligned_BSep_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag);
        
        %B01_raster_psth_sceneRecog_SAligned_BSep_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag);
        %B01_raster_psth_sceneRecog_SAligned_oldNew_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag);
        
        %B01_raster_psth_timeDiscrim_SAligned_BSep_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag);
        %B01_raster_psth_timeDiscrim_RAligned_BSep_NWB(  nwbData, timestampsOfCell, channelid, cellNr, brainAreaOfCell, bin_width_raster, bin_width_fr, plot_flag);
    end
end
# Code for paper Zheng et al. 'Hippocampal Theta Phase Precession Supports Memory Formation and Retrieval of Naturalistic Experience in Humans'

This repository contains the code for the following paper:

**Hippocampal Theta Phase Precession Supports Memory Formation and Retrieval of Naturalistic Experience in Humans**. Jie Zheng, Mar Yebra, Andrea G.P. Schjetnan, Clayton Mosher, Suneil K. Kalia, Jeffrey M. Chung, Chrystal M. Reed, Taufik A. Valiante, Adam N. Mamelak, Gabriel Kreiman, Ueli Rutishauser. *Nature Human Behavior, in press (2024)*. Pre-print (https://www.biorxiv.org/content/10.1101/2023.06.05.543539v1).

## Data
Neural recordings and behavior can be downloaded here (https://dandiarchive.org/dandiset/000940)

## Setup of code and Examples
** See (https://github.com/rutishauserlab/cogboundary-zheng) for how to setup the code and process/analyze the single-neuron data.

B01_raster_psth_encoding_BAligned_BSep_NWB: generate raster and PSTH plots for single neuron data during encoding, aligned to cognitive boundaries or virtual boundaries

B01_raster_psth_encoding_SAligned_BSep_NWB: generate raster and PSTH plots for single neuron data during encoding, aligned to clip onsets, separate for different boundary types

B01_raster_psth_sceneRecog_SAligned_BSep_NWB: generate raster and PSTH plots for single neuron data during scene recognition, aligned to image onsets, separate for different boundary types

B01_raster_psth_sceneRecog_SAligned_oldnew_NWB: generate raster and PSTH plots for single neuron data during scene recognition, aligned to image onsets, separate for target versus foil images 

B01_raster_psth_timeDiscrim_RAligned_BSep_NWB: genrate raster and PSTH plots for single neuron data during time discrimination, aligned to when participants respond, separate for different boundary types

B01_raster_psth_timeDiscrim_SAligned_BSep_NWB: generate raster and PSTH plots for single neuron data during time discrimination, aligned to image onsets, separate for different boundary types


** Codes for phase precession computation

C01_thetaBusrts_detection_NWB: detect theta bursts using the cycle-by-cycle analyses suite documented in Cole S and Voytek B, 2019, Journal of Neurophysiology (https://journals.physiology.org/doi/full/10.1152/jn.00273.2019)

C02_phase_precession_r_NWB: quantify phase precession by computing correlation coefficient between the spiking phases and time in unwrapped theta phases

C03_phase_precession_MI_NWB: quantify phase precession by computing the mutual information of spike-time in unwrapped theta phases relationship 

## Video clips
The accompanying video clips used in this task can be downloaded here (https://ucdavis.box.com/s/my3oitl4v35m195a19sjtswa34ozk9vd)


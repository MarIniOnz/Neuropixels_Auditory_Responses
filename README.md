# Neuropixels Auditory Responses
Lab Rotation Project performing Statistical Analysis on the responses of single and multi-units to auditory stimuli using Neuropixels 1.0 probe data, preprocessed with CatGT software (SpikeGLX feature).

## MATLAB codes.
 - **find_Clusters_Evoked_Response** : Get the spike times, clusters (single-units), and triggers, finding suitable clusters and plotting their frequency tuning, as long as the evoked responses for 8-tones trials in repeat.
 - **plot_Tuning_Curves** : Get the spike times, clusters and triggers. Plot clusters' tuning curves to different stimuli.
 - **view_Histogram_Dephts** : Get the spike times, clusters and triggers. PLotting drift map of the selected sequence, and compute histograms with the number of spikes that are contained within the specified depths.
 - **visualize_Plots** : Get the spike times, clusters and triggers. Choose a cluster and plot the raster plot of the spikes per trial.

## MATLAB Auxiliary codes.
  Auxiliary functions needed to run the main MATLAB codes.
  
## Python codes.
 - **NP_analysis_spikedata_with_triggers** : Get trigger times, spikes per trial and to which cluster they pertain.
 - **Spike_GLX_Analog_Trigger_loading_function** : Reading Analog triggers from AP bin file and from Analog bin file and convert them into a readable, operating variable.
 - **getTrigsMartin** : Reading Analog triggers from AP bin file and from Analog bin file and convert them into a readable, operating variable. Modified version of the previous one, more automated one, needing less inputs.

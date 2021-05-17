Analysis
------------------------------
Demixing, detrending, hemodynamics removal and parcellation - done on farnam
Input to this code: deltaF/F of parcels extracted by: Allen, Gal and a 4X4 grid
------------------------------
A. Trial data
1. Extract trials by main_time_traces_to_trials.m. Outputs are trial imaging and behavioral data
saved as: X:\Hadas\Meso-imaging\lan\xwpsych\spt\xw16imaging_time_traces_global.mat

2. Extract arousal state per trial and concatenate trials by state: pre_process_trials_animals.m
   Outputs are saved as:
X:\Hadas\Meso-imaging\lan\xwpsych\spt\xw_trial_imaging_time_traces_global_Allen.mat,
xw_trial_imaging_time_traces_global_LSSC, xw_trial_imaging_time_traces_global_Grid4,  and  xw_trial_meta_data.

  

3. Behavior Prediction by slope - main_behavior_prediction.m outputs are saved in 
   X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\behavior_prediction\xw_low_pup_q.mat
 
4. Network analysis by state per animal main_network_centrality_evaluateion_trials.m
   Output saved as: X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality\xx_low_pup_qtrials_incorrect_Gal.mat



------------------
B. Spont data
1. Extract spont time traces by main_time_traces_to_spont.m
   Outputs saved as: X:\Hadas\Meso-imaging\lan\xzpsych\spt\xz18imaging_time_traces_global_ITI.mat
2. Extract 3 sustained arousal states. Outputs saves as: 
   time stamps of segments: X:\Hadas\Meso-imaging\lan\xxpsych\spike2Features\xx12arousal_state_ITI_segemts.mat
3. Extract time traces per animal per state (for all days) pre_process_spont_animals.m
   Output saved as:  X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\xx_spont_data_3states.mat  
4. Network analysis by state per animal main_network_centrality_evaluateion_spont.m
   Output saved as: X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality\xs_partial_corr_high_pup_lAllen.mat

------------------
C. Continuous network analysis
1. main_network_to_arousal_LAN.m

Make figs - 

Figure 1 - behavior            - main_make_behavior_figs.m
Figure 2 - activity            - main_make_activity_figs.m
Figure 3 - behavior prediction - main_make_prediction_figs.m
Figure 4 - Network             - main_make_network_figs.m




Analysis


A. Trial data
1. Extract arousal state per trial and concatenate trials by state: pre_process_trials_animals_2p_lan.m
   Outputs are saved as:
time traces of imaging at:
X:\Hadas\Meso-imaging\lan\zapsych\spt\za_trial_imaging_time_traces_global_dfff.
trials labels (contrast, arousal, behavior, wheel, pupil):
X:\Hadas\Meso-imaging\lan\zapsych\spt\za_trial_meta_data.mat

3. Behavior Prediction by slope - main_behavior_prediction_2p_lan.m outputs are saved in 
   X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\behavior_prediction\za_low_pup_q.mat
4. Network analysis by state per animal main_network_centrality_evaluateion_trials_2p_lan.m
   Output saved as: X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\network_centrality\xx_low_pup_qtrials_incorrect_Gal.mat



------------------
B. Spont data
1. Extract spont time traces by main_time_traces_to_spont_2p_lan.m
   Outputs saved as: X:\Hadas\Meso-imaging\lan\xzpsych\spt\xz18imaging_time_traces_global_ITI.mat
2. Extract 3 sustained arousal states. Outputs saves as: 
   time stamps of segments: X:\Hadas\Meso-imaging\lan\xxpsych\spike2Features\xx12arousal_state_ITI_segemts.mat
3. Extract time traces per animal per state (for all days) pre_process_spont_animals_2p_lan.m
   Output saved as:  X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\xx_spont_data_3states.mat  
4. Network analysis by state per animal main_network_centrality_evaluateion_spont_2p_lan.m
   Output saved as: X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\network_centrality\xs_partial_corr_high_pup_lAllen.mat

------------------

Make figs - 

Figure 1 - behavior            - main_make_behavior_figs_2p_lan.m
Figure 2 - activity            - main_make_activity_figs_2p_lan.m
Figure 3 - behavior prediction - main_make_prediction_figs_2p_lan.m
Figure 4 - Network             - main_make_network_figs_2p_lan.m




to make a meso m file into a cripsr mfile: 
saving as _crispt
save results+figs into _crispr
delete all trial stuff

Analysis - Spont data
------------------------------
Demixing, detrending, hemodynamics removal and parcellation - done on farnam
Input to this code: deltaF/F of parcels extracted by: Allen and a 4X4 grid
------------------------------
1. Extract sustained states from smrx files: pre_process_spont_animals_crispr.m. 
   Outputs saved to data folder as arousal_state_ITI_segemts.mat

2. Extract time traces per animal per state (for all days) pre_process_spont_animals_crispr.m
   Output saved as:  X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr\...\arousal_state_ITI_segemts.mat  
3. Network analysis by state per animal main_network_centrality_evaluateion_spont_crispr.m
   Output saved as: X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality\xs_partial_corr_high_pup_lAllen.mat

------------------




Make figs - saved in X:\Hadas\Meso-imaging\crispr\meso_results\figs\
Figure 1 - behavior - main_make_behavior_figs_crispr.m
                      Uses arousal_state_ITI_segemts files (step 1 of analysis)

!!Figure 2 - activity            - main_make_activity_figs.m
!!Figure 4 - Network             - main_make_network_figs.m








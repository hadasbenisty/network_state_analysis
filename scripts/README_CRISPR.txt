to make a meso m file into a cripsr mfile: 
saving as _crispt
save results+figs into _crispr
delete all trial stuff

Analysis - Spont data
------------------------------
Demixing, detrending, hemodynamics removal and parcellation - done on farnam
Input to this code: deltaF/F of parcels extracted by: Allen and a 4X4 grid
------------------------------
0. Process smrx files: spike2processing_crispr.m results are save in 
   X:\Hadas\Meso-imaging\Antara\data\Antara\AnalyzedData\...\smrx_signals_v4.mat
1. Extract sustained states from smrx files: pre_process_spont_animals_crispr.m. 
   Outputs saved to data folder as arousal_traces_states.mat
   Output saved as:  X:\Hadas\Meso-imaging\CRISPR\analysis_results\...\con_states.mat
2. Correlation matrices & network analysis by state per animal main_network_centrality_evaluateion_spont_crispr.m
   Output saved as: X:\Hadas\Meso-imaging\CRISPR\analysis_results\network_centrality_pearson_corr\

------------------




Make figs - saved in X:\Hadas\Meso-imaging\crispr\meso_results\figs\
Figure 0 - genral report per animal (spont session) - main_plot_initial_report.m
Figure 1 - arousal                 - main_make_arousal_figs_crispr.m                      
Figure 2 - activity                 - main_make_activity_figs_crispr.m
Figure 3 - Correlation & Network    - main_make_network_figs_crispr.m








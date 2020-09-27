%%Code to process the ITI period for all animals with extraction 
%of running periods, quiescence periods, and high pupil and low pupil periods.
restoredefaultpath;
rehash toolboxcache
cd('X:\Hadas\Meso-imaging\lan\results\code\Directed_Network_Allen')
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions'));
main_organize_by_ITI_global_parcellation_lan({'xx'}, [12 14 15 16 17 18 19 20 21 22 23])
main_organize_by_ITI_global_parcellation_lan({'xx'}, [25 26 27 29 30 33 34 35])
main_organize_by_ITI_global_parcellation_lan({'xw'}, [13 14 15 16 17 18 19 21 22 23])
main_organize_by_ITI_global_parcellation_lan({'xw'}, [24 25 26 27 28 29 30 33 34 35 36])
main_organize_by_ITI_global_parcellation_lan({'xz'}, [18 19 20 21 22 23 24 25])
main_organize_by_ITI_global_parcellation_lan({'xz'}, [26 27 28 29 30 33 34 35 36])
main_organize_by_ITI_global_parcellation_lan({'xu'}, [20 22 24 26 27 29])
main_organize_by_ITI_global_parcellation_lan({'xu'}, [30 32 33 34 35 36 37])
main_organize_by_ITI_global_parcellation_lan({'xs'}, [15 16 17 18 19 20])
main_organize_by_ITI_global_parcellation_lan({'xs'}, [21 24 25 27 29 30 31 32])
main_organize_by_ITI_global_parcellation_lan({'xt'}, [15 17 18 19 21 22 23])
main_organize_by_ITI_global_parcellation_lan({'xt'}, [25 26 27 28 29 30 32])

%check

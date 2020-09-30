%% spon high low pupil
cd 'D:\network_state_analysis\scripts\parcor_undirected\spon_highpup_lowpup';
step1_weightsandgraphanalysis
cd 'D:\network_state_analysis\scripts\parcor_undirected\spon_highpup_lowpup';
outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory_spet29\';
step2_averageplots(inputpth,outputpth);
cd 'D:\network_state_analysis\scripts\parcor_undirected\spon_highpup_lowpup';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots(inputpth,outputpth);

%% spon run not
cd 'D:\network_state_analysis\scripts\parcor_undirected\spon_run_notrun';
step1_weightsandgraphanalysis
cd 'D:\network_state_analysis\scripts\parcor_undirected\spon_run_notrun';
outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory_spet29\';
step2_averageplots(inputpth,outputpth);

cd 'D:\network_state_analysis\scripts\parcor_undirected\spon_run_notrun';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots(inputpth,outputpth);


%% trial
cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
step1_weightsandgraphanalysis

%% 1
cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory_spet29\';
step2_averageplots_lowpupil_corrincorr(inputpth,outputpth)

cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_lowpupil_corrincorr(inputpth,outputpth)

%% 2
cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory_spet29\';
step2_averageplots_notrun_corrincorr(inputpth,outputpth)

cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_notrun_corrincorr(inputpth, outputpth)

%% 3
cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory_spet29\';
step2_averageplots_run_corrincorr(inputpth,outputpth)

cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_run_corrincorr(inputpth, outputpth)

%% 4
cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory_spet29\';
step2_averageplots_highpupil_corrincorr(inputpth,outputpth)

cd 'D:\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_highpupil_corrincorr(inputpth, outputpth)





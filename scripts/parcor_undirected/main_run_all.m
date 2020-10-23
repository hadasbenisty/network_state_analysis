addpath(genpath('X:\Lav\network_state_analysis\utils'))
%% spon high low pupil
cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_highpup_lowpup';
step1_weightsandgraphanalysis
% cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_highpup_lowpup';
% outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
% inputpth='X:\Lav\ProcessingDirectory_spet29\';
% step2_averageplots(inputpth,outputpth);
cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_highpup_lowpup';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots(inputpth,outputpth);
%% spon high low face
%cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_highface_lowface';
%step1_weightsandgraphanalysis
%cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_highface_lowface';
%outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
%inputpth='X:\Lav\ProcessingDirectory\';
%step2_averageplots(inputpth,outputpth);

%% spon run not
cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_run_notrun';
step1_weightsandgraphanalysis
% cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_run_notrun';
%outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
% inputpth='X:\Lav\ProcessingDirectory_spet29\';
% step2_averageplots(inputpth,outputpth);

cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_run_notrun';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots(inputpth,outputpth);


%% trial
% cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
% step1_weightsandgraphanalysis

% %% 1
% cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
% outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
% inputpth='X:\Lav\ProcessingDirectory_spet29\';
% step2_averageplots_lowpupil_corrincorr(inputpth,outputpth)

cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_lowpupil_corrincorr(inputpth,outputpth)

% %% 2
% cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
% outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
% inputpth='X:\Lav\ProcessingDirectory_spet29\';
% step2_averageplots_notrun_corrincorr(inputpth,outputpth)

cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_notrun_corrincorr(inputpth, outputpth)

%% 3
% cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
% outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
% inputpth='X:\Lav\ProcessingDirectory_spet29\';
% step2_averageplots_run_corrincorr(inputpth,outputpth)

cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_run_corrincorr(inputpth, outputpth)

% %% 4
% cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
% outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
% inputpth='X:\Lav\ProcessingDirectory_spet29\';
% step2_averageplots_highpupil_corrincorr(inputpth,outputpth)

cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\beforetrial_corr_incorr';
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_highpupil_corrincorr(inputpth, outputpth)

cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_3states'
outputpth = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
inputpth='X:\Lav\ProcessingDirectory\';
step2_averageplots_3states(inputpth, outputpth)

% cd 'X:\Lav\network_state_analysis\scripts\parcor_undirected\spon_3states'
% outputpth = 'X:\Lav\ProcessingDirectory_spet29\parcor_undirected\';
% inputpth='X:\Lav\ProcessingDirectory_spet29\';
% step2_averageplots_3states(inputpth, outputpth)





function main_network_centrality_evaluateion_trials_2p_lan
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
addpath(genpath('../network_state_analysis\gspbox'));
T = readtable('X:\Lan\FMB208 conditioning imaging data\animal_log_imaging_lt.xlsx');
animals = T.AnimalID;
dayslist = T.PsychTest;
i = find(strcmp(animals, 'zb'));
i(2) = find(strcmp(animals, 'zy'));
i(3) = find(strcmp(animals, 'zn'));
i(4) = find(strcmp(animals, 'zi'));
i(5) = find(strcmp(animals, 'zc'));
i(6) = find(strcmp(animals, 'zm'));

animals=animals(setdiff(1:length(animals), i));
dayslist=dayslist(setdiff(1:length(dayslist), i));

pre_trial_time_start = -3;
pre_trial_time_end = -.1;
doover=1;maxdays=30;
similarity_name = { 'pearson_corr'  };%'corr',,  'fullcorr' 'cov','partial_corr'
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
for sim_i = 1:length(similarity_name)
    for ai = 1:length(animals)
        eval_weights_and_cent_per_day(eval(dayslist{ai}), maxdays, doover, similarity_name{sim_i}, animals{ai}, statenames, pre_trial_time_start, pre_trial_time_end, 0);
        
        %      eval_weights_and_cent(maxdays, doover, similarity_name{sim_i}, animals{ai}, statenames, pre_trial_time_start, pre_trial_time_end, 0);
        %     for permi = 1:100
        %         eval_weights_and_cent(isloose, animals{ai}, statenames, pre_trial_time_start, pre_trial_time_end, permi);
        %     end
    end
end
end
function [data, data_corr, data_inco, slopes_corr, slopes_inco, slopes_corr_baseline, slopes_inco_baseline] = ...
    get_trial_data(animal, statename, pre_trial_time_start, pre_trial_time_end, toperm, maxdays, daynum)
data = [];
data_corr = [];
data_inco = [];
slopes_corr=[]; slopes_inco=[];
slopes_corr_baseline=[]; slopes_inco_baseline=[];
if ~exist(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_meta_data.mat'], 'file')
    
    return;
end
load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_meta_data.mat'], ...
    'behavior_labels', 'arousal_labels', 'days_labels');
sel_trials = behavior_labels < 3 & arousal_labels == statename2num( statename)&days_labels<=maxdays;
if exist('daynum','var')
    sel_trials = sel_trials & days_labels == daynum;
end
suc_fail_labels = behavior_labels(sel_trials);
if toperm
    suc_fail_labels = suc_fail_labels(randperm(length(suc_fail_labels)));
end
load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_imaging_time_traces_global_dfff.mat'], 'imaging_time_tracesN', 't');

imaging_time_tracesN= imaging_time_tracesN(:, :, sel_trials);
imaging_time_traces_cor = imaging_time_tracesN(:, :, suc_fail_labels==1);
imaging_time_traces_inc = imaging_time_tracesN(:, :, suc_fail_labels==2);


data=[];data_inco=[];data_corr=[];
slopes_corr=[];slopes_inco=[];
for T=1:size(imaging_time_tracesN,3)
    if all(~isnan(imaging_time_tracesN(:,t>=pre_trial_time_start & t<pre_trial_time_end,T)))
        data = cat(2, data,  imaging_time_tracesN(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
end
for T=1:size(imaging_time_traces_cor,3)
    if all(all(~isnan(imaging_time_traces_cor(:,t>=pre_trial_time_start & t<pre_trial_time_end,T))))
        data_corr = cat(2, data_corr,  imaging_time_traces_cor(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
        slopes_corr(:,T) = getSlope(imaging_time_traces_cor(:, t>0 & t>0.33,T));
        slopes_corr_baseline(:,T) = getSlope(imaging_time_traces_cor(:, t>-0.33 & t>0,T));
    end
end
for T=1:size(imaging_time_traces_inc,3)
    if all(~isnan(imaging_time_traces_inc(:,t>=pre_trial_time_start & t<pre_trial_time_end,T)))
        data_inco = cat(2, data_inco,  imaging_time_traces_inc(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
        slopes_inco(:,T) = getSlope(imaging_time_traces_inc(:, t>0 & t>0.33,T));
        slopes_inco_baseline(:,T) = getSlope(imaging_time_traces_inc(:, t>-0.33 & t>0,T));
    end
end

end
function eval_weights_and_cent_per_day(dayslist, maxdays, doover, simname, animal, statenames, pre_trial_time_start, pre_trial_time_end, topermind)

outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\network_centrality_' simname];
mkNewDir(outputfolder);
signalnames = {'cells' };%'grid4' ,'LSSC''Allen'
disp(animal)
for sig_i = 1:length(signalnames)
    for day_i = 1:length(dayslist)
        for state_i = 1:length(statenames)
            disp([animal ' ' num2str(dayslist(day_i)) ' ' statenames{state_i}])
            [~, data_corr, data_inco, slopes_corr, slopes_incorr, slopes_corr_baseline, slopes_inco_baseline] = ...
                get_trial_data(animal, statenames{state_i}, pre_trial_time_start, pre_trial_time_end, topermind, maxdays, dayslist(day_i));
            ththird = round(size(data_corr,1));
            for th= ththird%3:2:51
                if exist('W_corr_cor','var')
                    if th > size(W_corr_cor,1)
                        return;
                    end
                end
                corrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} '_D' num2str(dayslist(day_i)), '_ththird.mat']);
                    incorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} '_D' num2str(dayslist(day_i)), '_ththird.mat']);
                    
%                 if topermind == 0
%                     corrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} '_D' num2str(dayslist(day_i)), '_' num2str(th) '.mat']);
%                     incorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} '_D' num2str(dayslist(day_i)), '_' num2str(th) '.mat']);
%                     
%                 else
%                     corrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} 'perm' num2str(topermind) '_D' num2str(dayslist(day_i)), '_' num2str(th) '.mat']);
%                     incorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} 'perm' num2str(topermind) '_D' num2str(dayslist(day_i)), '_' num2str(th) '.mat']);
%                 end
                if exist(corrfile, 'file') && exist(incorrfile, 'file')&&~doover
                    load(corrfile, 'W_corr_cor');
                    load(incorrfile,'W_corr_inc');
                else
                    if isempty(data_corr)
                        continue;
                    end
                    
                    W_corr_cor = measure_weights(data_corr, simname);
                    W_corr_inc = measure_weights(data_inco, simname);
                end
                
                % correct
                if ~isempty(W_corr_cor)
                    [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_cor, 1:size(W_corr_cor), [], @process_sim, th);
                    save(corrfile,'W_corr_cor',...
                        'cent_corr_weighted',...
                        'cent_corr_notweighted', 'G_corr', 'names_corr', 'slopes_corr_baseline','slopes_corr');
                end
                
                % incorrect
                if ~isempty(W_corr_inc)
                    [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_inc, 1:size(W_corr_inc), [], @process_sim, th);
                    save(incorrfile,'W_corr_inc',...
                        'cent_corr_weighted',...
                        'cent_corr_notweighted', 'G_corr', 'names_corr', 'slopes_inco_baseline','slopes_incorr');
                end
            end
        end
    end
end
end

function eval_weights_and_cent(maxdays, doover, simname, animal, statenames, pre_trial_time_start, pre_trial_time_end, topermind)

outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\network_centrality_' simname];
mkNewDir(outputfolder);
signalnames = {'cells' };%'grid4' ,'LSSC''Allen'
disp(animal)
for sig_i = 1:length(signalnames)
    for state_i = 1:length(statenames)
        disp(statenames{state_i})
        for th= 3:2:51
            if exist('W_corr_cor','var')
                if th > size(W_corr_cor,1)
                    return;
                end
            end
            if topermind == 0
                corrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} '_' num2str(th) '.mat']);
                incorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} '_' num2str(th) '.mat']);
                
            else
                corrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} 'perm' num2str(topermind) '_' num2str(th) '.mat']);
                incorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} 'perm' num2str(topermind) '_' num2str(th) '.mat']);
            end
            if exist(corrfile, 'file') && exist(incorrfile, 'file')&&~doover
                load(corrfile, 'W_corr_cor');
                load(incorrfile,'W_corr_inc');
            else
                [~, data_corr, data_inco] = ...
                    get_trial_data(animal, statenames{state_i}, pre_trial_time_start, pre_trial_time_end, topermind, maxdays);
                if isempty(data_corr)
                    continue;
                end
                W_corr_cor = measure_weights(data_corr, simname);
                W_corr_inc = measure_weights(data_inco, simname);
            end
            
            % correct
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_cor, 1:size(W_corr_cor), [], @process_sim, th);
            save(corrfile,'W_corr_cor',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
            % incorrect
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_inc, 1:size(W_corr_inc), [], @process_sim, th);
            save(incorrfile,'W_corr_inc',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
        end
    end
end
end


% function [corralldata,incorralldata]=permute_corrincorr_lan(imaging_time_traces_cor,imaging_time_traces_inc)
% %concatenate conditions
% concatenated_dff=cat(3,imaging_time_traces_cor,imaging_time_traces_inc);
% %size of each condition and size of the entire dataset to permute
% samplesize=size(imaging_time_traces_cor,3);
% permutation = randperm(size(concatenated_dff,3));
% %save permutations
% corralldata=concatenated_dff(:,:,permutation(1:samplesize));
% incorralldata=concatenated_dff(:,:,permutation(samplesize+1:length(permutation)));
% %save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',num2str(i),animal,'perm_running_time_traces'),'runningalldata','runningalldata_t','notrunningalldata_t','notrunningalldata');
% end
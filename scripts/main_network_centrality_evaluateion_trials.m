function main_network_centrality_evaluateion_trials
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
addpath(genpath('../network_state_analysis\gspbox'));
animals={'xt','xu' 'xs', 'xx','xz','xw'};%,
pre_trial_time_start = -5;
pre_trial_time_end = -.1;
doover=1;maxdays=30;
similarity_name = { 'pearson_corr'};
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
for sim_i = 1:length(similarity_name)
for ai = 1:length(animals)
%    eval_weights_and_cent_per_day(animals{ai}, statenames, pre_trial_time_start, pre_trial_time_end, 0);

     eval_weights_and_cent(maxdays, doover, similarity_name{sim_i}, animals{ai}, statenames, pre_trial_time_start, pre_trial_time_end, 0);
    %     for permi = 1:100
    %         eval_weights_and_cent(isloose, animals{ai}, statenames, pre_trial_time_start, pre_trial_time_end, permi);
    %     end
end
end
end
function [data, data_corr, data_inco, roi_names, region_labels] = ...
    get_trial_data(signame, animal, statename, pre_trial_time_start, pre_trial_time_end, toperm, maxdays, daynum)
roi_names=[];
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
load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_imaging_time_traces_global_' signame '_dfff.mat'], 'imaging_time_traces', 't',...
    'roi_names','region_labels');

imaging_time_traces= imaging_time_traces(:, :, sel_trials);
imaging_time_traces_cor = imaging_time_traces(:, :, suc_fail_labels==1);
imaging_time_traces_inc = imaging_time_traces(:, :, suc_fail_labels==2);


data=[];data_inco=[];data_corr=[];

for T=1:size(imaging_time_traces,3)
    if all(~isnan(imaging_time_traces(:,t>=pre_trial_time_start & t<pre_trial_time_end,T)))
        data = cat(2, data,  imaging_time_traces(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
end
for T=1:size(imaging_time_traces_cor,3)
    if all(~isnan(imaging_time_traces_cor(:,t>=pre_trial_time_start & t<pre_trial_time_end,T)))
        data_corr = cat(2, data_corr,  imaging_time_traces_cor(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
end
for T=1:size(imaging_time_traces_inc,3)
    if all(~isnan(imaging_time_traces_inc(:,t>=pre_trial_time_start & t<pre_trial_time_end,T)))
        data_inco = cat(2, data_inco,  imaging_time_traces_inc(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
end

end
function eval_weights_and_cent_per_day(animal, statenames, pre_trial_time_start, pre_trial_time_end, topermind)

outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory';
% [~, allen_parcels] = getParcellsByLansAllansAtlas;
% [parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
%
% regionLabel.Allen = allen_parcels.regionNum;
% regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
% [parcels_names.Gal, ~, ~, regionLabel.Gal] = get_gal_meta_parcels_by_allen(parcels_names.Allen, maskByAllen.Allen, ...
%     regionLabel.Allen, animal);
signalnames = {'Allen','LSSC'};
disp(animal)
[~,days_list] = animaltodays(animal);
for day_i = 1:length(days_list)
    for sig_i = 1:length(signalnames)
        for state_i = 1:length(statenames)
            disp(statenames{state_i})
            disp([animal ' ' statenames{state_i}  ' ' signalnames{sig_i} ' day ' num2str(days_list(day_i)) ]);
            if topermind == 0
                corrfile = fullfile(outputfolder,'network_centrality',[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} '_day' num2str(days_list(day_i)) '.mat']);
                incorrfile = fullfile(outputfolder,'network_centrality',[animal '_' statenames{state_i} ,'trials_incorrect_' signalnames{sig_i} '_day' num2str(days_list(day_i)) '.mat']);
                
            else
                corrfile = fullfile(outputfolder,'network_centrality',[animal '_' statenames{state_i} ,'trials_correct_' signalnames{sig_i} 'perm' num2str(topermind) '_day' num2str(days_list(day_i)) '.mat']);
                incorrfile = fullfile(outputfolder,'network_centrality',[animal '_' statenames{state_i} ,'trials_incorrect_' signalnames{sig_i} 'perm' num2str(topermind) '_day' num2str(days_list(day_i)) '.mat']);
            end
            if exist(corrfile, 'file') && exist(incorrfile, 'file')
                load(corrfile, 'W_corr_cor');
                load(incorrfile,'W_corr_inc');
                [~, ~, ~, parcels_names, regionLabel] = ...
                    get_trial_data(signalnames{sig_i}, animal, statenames{state_i}, pre_trial_time_start, pre_trial_time_end, topermind);
                
            else
                [~, data_corr, data_inco, parcels_names, regionLabel] = ...
                    get_trial_data(signalnames{sig_i}, animal, statenames{state_i}, pre_trial_time_start, pre_trial_time_end, topermind, days_list(day_i));
                W_corr_cor = measure_weights_partial(data_corr, 'corr');
                W_corr_inc = measure_weights_partial(data_inco, 'corr');
            end
            % correct
            if all(~isnan(W_corr_cor(:)))&& ~isempty(W_corr_cor)
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_cor, parcels_names, regionLabel);
            save(corrfile,'W_corr_cor',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
            end
            % incorrect
             if all(~isnan(W_corr_inc(:)))&& ~isempty(W_corr_inc)
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_inc, parcels_names, regionLabel);
            save(incorrfile,'W_corr_inc',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
             end
        end
    end
end
end
function eval_weights_and_cent(maxdays, doover, simname, animal, statenames, pre_trial_time_start, pre_trial_time_end, topermind)

outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality_' simname];
mkNewDir(outputfolder);
signalnames = {'Allen'};% 'LSSC' 'Grid4'  };%'Allen' 
disp(animal)
for sig_i = 1:length(signalnames)
    switch signalnames{sig_i} 
        case {'Allen','LSSC'}
            thvals =[inf 7];%5:2:23];
        case 'Grid4'
            thvals = [inf 150];% 150:50:400];
    end
    for state_i = 1:length(statenames)
        disp(statenames{state_i})
        for th= thvals
        if topermind == 0
            Wcorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} '_W.mat']);
            Wincorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} '_W.mat']);
            corrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} '_' num2str(th) '.mat']);
            incorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} '_' num2str(th) '.mat']);
            
        else
            Wcorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} 'perm' num2str(topermind) '_W.mat']);
            Wincorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} 'perm' num2str(topermind) '_W.mat']);
            corrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_correct_' signalnames{sig_i} '_' num2str(th) '.mat']);
            incorrfile = fullfile(outputfolder,[animal '_' statenames{state_i} ,'_trials_incorrect_' signalnames{sig_i} '_' num2str(th) '.mat']);
            
        end
        if exist(Wcorrfile, 'file') && exist(Wincorrfile, 'file')&&~doover
            load(Wcorrfile, 'W_corr_cor');
            load(Wincorrfile,'W_corr_inc');
            roi_names=[];
            load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_imaging_time_traces_global_' signalnames{sig_i} '_dfff.mat'], ...
                'roi_names','region_labels');
            parcels_names=roi_names;
        else
            [~, data_corr, data_inco, parcels_names, region_labels] = ...
                get_trial_data(signalnames{sig_i}, animal, statenames{state_i}, pre_trial_time_start, pre_trial_time_end, topermind, maxdays);
            W_corr_cor = measure_weights(data_corr, simname);
            W_corr_inc = measure_weights(data_inco, simname);
            save(Wcorrfile, 'W_corr_cor');
            save(Wincorrfile,'W_corr_inc');
        end
        
            % correct
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_cor, parcels_names, region_labels, @process_sim, th);
            save(corrfile,'W_corr_cor',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
            % incorrect
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_inc, parcels_names, region_labels, @process_sim, th);
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
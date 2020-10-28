function main_network_centrality_evaluateion_trials
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
saveplots = false;
animals={'xs','xx','xz','xw','xt','xu'};
pre_trial_time_start = -3;
pre_trial_time_end = -.1;
isloose = true;
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
% for ai = 1:length(animals)
%     eval_weights_and_cent(isloose, animals{ai}, saveplots, statenames, pre_trial_time_start, pre_trial_time_end);
% end
for ismidcontrast=0:1
outputfiggolder = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
% [trials_states_notweighted, trials_states_weighted] = plot_centrality_res(ismidcontrast, isloose, animals, outputfiggolder, statenames, 'trials');
% [correct_states_notweighted, correct_states_weighted] = plot_centrality_res(ismidcontrast, isloose, animals, outputfiggolder, statenames, 'trials_correct');
% [incorrect_states_notweighted, incorrect_states_weighted] = plot_centrality_res(ismidcontrast, isloose, animals, outputfiggolder, statenames, 'trials_incorrect');
close all;
if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
if ismidcontrast
    midcontraststr = 'mid_contrast';
else
    midcontraststr='';
end
% save(fullfile(outputfiggolder, ['centrality_stats_pretrial' loosestr midcontraststr '.mat']), 'trials_states_notweighted',...
%     'trials_states_weighted', 'correct_states_notweighted', 'correct_states_weighted',...
%     'incorrect_states_notweighted', 'incorrect_states_weighted');
%plotSummaryCentrality(ismidcontrast, isloose, outputfiggolder, statenames);
makeslopeamplitudeplots(animals, isloose,outputfiggolder)
end
end
%%
function plotSummaryCentrality(ismidcontrast, isloose, outputfiggolder, statenames)
if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
if ismidcontrast
    midcontraststr = 'mid_contrast';
else
    midcontraststr='';
end
load(fullfile(outputfiggolder, ['centrality_stats_pretrial' loosestr  midcontraststr '.mat']), ...
    'correct_states_notweighted', 'correct_states_weighted',...
    'incorrect_states_notweighted', 'incorrect_states_weighted');
[parcels_names] = get_allen_meta_parcels;
braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
N = size(correct_states_notweighted.low_pup_q.eigenvector,2);
%% not weighted
centnames = fieldnames(correct_states_notweighted.low_pup_q);
for l=1:length(centnames)
    centstr = centnames{l};
    figure;
for k=1:length(statenames) 
M1 = nanmean(correct_states_notweighted.(statenames{k}).(centstr),2);
S1 = nanstd(correct_states_notweighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
M2 = nanmean(incorrect_states_notweighted.(statenames{k}).(centstr),2);
S2 = nanstd(incorrect_states_notweighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
subplot(3, 1, k);
barwitherr([S1 S2], [M1 M2]);
set(gca,'XTick', 1:length(parcels_names))
set(gca,'XTickLabel', parcels_names)
strttl = statenames{k};
strttl(strttl=='_') = ' ';
title(strttl);
axis tight;
end
legend('Correct','Incorrect');
suptitle(['Correct/Incorrect ' centstr]);
set(gcf, 'Position', [1          41        1920         963]);
mysave(gcf, fullfile(outputfiggolder, 'not_weighted', [centnames{l} '_centrality_stats_pretrial' loosestr midcontraststr]));
end
%%  weighted
centnames = fieldnames(correct_states_weighted.low_pup_q);
for l=1:length(centnames)
    centstr = centnames{l};
    figure;
for k=1:length(statenames) 
M1 = nanmean(correct_states_weighted.(statenames{k}).(centstr),2);
S1 = nanstd(correct_states_weighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
M2 = nanmean(incorrect_states_weighted.(statenames{k}).(centstr),2);
S2 = nanstd(incorrect_states_weighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
subplot(3, 1, k);
barwitherr([S1 S2], [M1 M2]);
set(gca,'XTick', 1:length(parcels_names))
set(gca,'XTickLabel', parcels_names)
strttl = statenames{k};
strttl(strttl=='_') = ' ';
title(strttl);
axis tight;
end
legend('Correct','Incorrect');
suptitle(['Correct/Incorrect ' centstr]);
set(gcf, 'Position', [1          41        1920         963]);
mysave(gcf, fullfile(outputfiggolder, 'weighted', [centnames{l} '_centrality_stats_pretrial' loosestr  midcontraststr]));
%% difference
centnames = fieldnames(correct_states_notweighted.low_pup_q);
for centt=1:length(centnames)
    centstr = centnames{centt};
    for k=1:length(statenames)
        M1 = nanmean(correct_states_notweighted.(statenames{k}).(centstr),2);
        M2 = nanmean(incorrect_states_notweighted.(statenames{k}).(centstr),2);
        
        %make graphs
        graph_overlay_allen_paired([loosestr midcontraststr],fullfile(outputfiggolder, 'not_weighted'), M1,...
            M2,'trials',strcat(statenames{k},'_corrincorr_',centstr),['corr - incorr ' centstr ' Centrality (trials)' statenames{k}],parcels_names,N);
        
        graph_heatmap([loosestr midcontraststr],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
            M1,M2,'trials',strcat(statenames{k},'_corrincorr_',centstr),['corr - incorr ' centstr ' Centrality (trials)' statenames{k}]);
    end
end
end
end
%%
function [spon_states_notweighted, spon_states_weighted] = plot_centrality_res(ismidcontrast, isloose, animals, outputfiggolder, statenames, suffix_files)
if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
if ismidcontrast
    midcontraststr = 'mid_contrast';
else
    midcontraststr='';
end
[parcels_names] = get_allen_meta_parcels;
cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector'};
for state_i = 1:length(statenames)
    for cent_i = 1:length(cent_features)
        spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = [];
        spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = [];
    end
end
for i=1:length(animals)
    animal=animals{i};
    outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020',animal);
    for state_i = 1:length(statenames)
        
        
        load(fullfile(outputfolder,['network_analysis_corr',statenames{state_i} ,suffix_files, loosestr, midcontraststr, '.mat']), ...
            'cent_corr_weighted','cent_corr_notweighted');
        cent_features = fieldnames(cent_corr_weighted);
        for cent_i = 1:length(cent_features)
            spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
                spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}), ...
                cent_corr_notweighted.(cent_features{cent_i}));
            
            spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
                spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}), ...
                cent_corr_weighted.(cent_features{cent_i}));
        end
        
    end
end
mkNewDir(fullfile(outputfiggolder, 'weighted'))
mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
legstr = {'Low Q', 'High Q', 'Loc'};
for ni = 1:length(cent_features)
  
    graph_overlay_allen_3conditions([loosestr midcontraststr], fullfile(outputfiggolder, 'weighted'), spon_states_weighted.low_pup_q.(cent_features{ni}),...
        spon_states_weighted.high_pup_q.(cent_features{ni}), spon_states_weighted.high_pup_l.(cent_features{ni}),...
        suffix_files,cent_features{ni},['3 states ' cent_features{ni} ' Centrality (Trials)'], parcels_names,length(animals), legstr);
    
    graph_overlay_allen_3conditions([loosestr midcontraststr], fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.low_pup_q.(cent_features{ni}),...
        spon_states_notweighted.high_pup_q.(cent_features{ni}), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
        suffix_files,cent_features{ni},['3 states ' cent_features{ni} ' Centrality (Trials)'], parcels_names,length(animals), legstr);
    
    %states as before
    mkNewDir(fullfile(outputfiggolder, 'not_weighted',suffix_files,'spon_pupilhigh_pupillow'))
    mkNewDir(fullfile(outputfiggolder, 'not_weighted',suffix_files,'spon_run_pupillow'))
    mkNewDir(fullfile(outputfiggolder, 'not_weighted',suffix_files,'spon_run_pupilhigh'))
    
    %difference maps, not weighted, for each centrality measure
    braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
    parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');

    %high vs low pup
    graph_overlay_allen_paired([loosestr midcontraststr],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_q.(cent_features{ni}),...
        spon_states_notweighted.low_pup_q.(cent_features{ni}),strcat(suffix_files,'/spon_pupilhigh_pupillow/'),cent_features{ni},['pupil high vs low ' cent_features{ni} ' Centrality (trials)'],parcels_names,length(animals));
    
    graph_heatmap([loosestr midcontraststr],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
        spon_states_notweighted.high_pup_q.(cent_features{ni}),spon_states_notweighted.low_pup_q.(cent_features{ni}),...
        strcat(suffix_files,'/spon_pupilhigh_pupillow/'),cent_features{ni},['pupil high vs low ' cent_features{ni} ' Centrality (trials)']);
    
    %locomotion vs low pup
    graph_overlay_allen_paired([loosestr midcontraststr],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
        spon_states_notweighted.low_pup_q.(cent_features{ni}),strcat(suffix_files,'/spon_run_pupillow/'),cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (trials)'],parcels_names,length(animals));

    graph_heatmap([loosestr midcontraststr],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
        spon_states_notweighted.high_pup_l.(cent_features{ni}),spon_states_notweighted.low_pup_q.(cent_features{ni}),...
        strcat(suffix_files,'/spon_run_pupillow/'),cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (trials)']);
    
    %locomotion vs high pup
    graph_overlay_allen_paired([loosestr midcontraststr],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
        spon_states_notweighted.high_pup_q.(cent_features{ni}),strcat(suffix_files,'/spon_run_pupilhigh/'),cent_features{ni},['run vs pupil high ' cent_features{ni} ' Centrality (trials)'],parcels_names,length(animals));
    
    graph_heatmap([loosestr midcontraststr],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
        spon_states_notweighted.high_pup_l.(cent_features{ni}),spon_states_notweighted.high_pup_q.(cent_features{ni}),...
        strcat(suffix_files,'/spon_run_pupilhigh/'),cent_features{ni},['run vs pupil high ' cent_features{ni} ' Centrality (trials)']);
 
end

end

function eval_weights_and_cent(isloose, animal, saveplots, statenames, pre_trial_time_start, pre_trial_time_end)
if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
mkNewDir(outputfolder);


load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states' loosestr '.mat'],...
    'low_pup_q','high_pup_q','high_pup_l','days_to_process', 't'); %#ok<NASGU>
[parcels_names, parcels_region_labels, finalindex] = get_allen_meta_parcels;
disp(animal)
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    if 0&&exist(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials.mat'),'file')
        continue;
    end
    data_3D = eval(statenames{state_i});
    imaging_time_traces_all = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary<3);
    imaging_time_traces_cor = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary==1);
    imaging_time_traces_inc = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary==2);
    imaging_time_traces_all_mid_contrast = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary<3&data_3D.trialslabels.contrastLabels>=10 & data_3D.trialslabels.contrastLabels <=40);
    imaging_time_traces_cor_mid_contrast = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary==1&data_3D.trialslabels.contrastLabels>=10 & data_3D.trialslabels.contrastLabels <=40);
    imaging_time_traces_inc_mid_contrast = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary==2&data_3D.trialslabels.contrastLabels>=10 & data_3D.trialslabels.contrastLabels <=40);
    
    
    data=[];data_inco=[];data_corr=[];
    data_mid_contrast=[];data_inco_mid_contrast=[];data_corr_mid_contrast=[];
    for T=1:size(imaging_time_traces_all,3)
       data = cat(2, data,  imaging_time_traces_all(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_cor,3)
       data_corr = cat(2, data_corr,  imaging_time_traces_cor(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_inc,3)
       data_inco = cat(2, data_inco,  imaging_time_traces_inc(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_all_mid_contrast,3)
       data_mid_contrast = cat(2, data_mid_contrast,  imaging_time_traces_all_mid_contrast(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_cor_mid_contrast,3)
       data_corr_mid_contrast = cat(2, data_corr_mid_contrast,  imaging_time_traces_cor_mid_contrast(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_inc_mid_contrast,3)
       data_inco_mid_contrast = cat(2, data_inco_mid_contrast,  imaging_time_traces_inc_mid_contrast(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
       
    data = data(finalindex, :);
    data_corr = data_corr(finalindex, :);
    data_inco = data_inco(finalindex, :);
    
    data_mid_contrast = data_mid_contrast(finalindex, :);
    data_corr_mid_contrast = data_corr_mid_contrast(finalindex, :);
    data_inco_mid_contrast = data_inco_mid_contrast(finalindex, :);
    
    data = data(:, all(~isnan(data)));
    data_corr = data_corr(:, all(~isnan(data_corr)));
    data_inco = data_inco(:, all(~isnan(data_inco)));
    
    data_mid_contrast = data_mid_contrast(:, all(~isnan(data_mid_contrast)));
    data_corr_mid_contrast = data_corr_mid_contrast(:, all(~isnan(data_corr_mid_contrast)));
    data_inco_mid_contrast = data_inco_mid_contrast(:, all(~isnan(data_inco_mid_contrast)));
    if any(isnan(data(:)))
        disp('nans in dataset')
        continue;
    end
    %% measure weights
    W_corr = measure_weights_partial(data, 'corr');
    W_corr_cor = measure_weights_partial(data_corr, 'corr');
    W_corr_inc = measure_weights_partial(data_inco, 'corr');
    
    W_corr_mid_contrast = measure_weights_partial(data_mid_contrast, 'corr');
    W_corr_cor_mid_contrast = measure_weights_partial(data_corr_mid_contrast, 'corr');
    W_corr_inc_mid_contrast = measure_weights_partial(data_inco_mid_contrast, 'corr');
    
    disp('corweights done')
    % Graph Analysis
    [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names);
   
    save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials', loosestr, '.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    
     [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_cor, parcels_names);
     save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials_correct', loosestr, '.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    
   [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_inc, parcels_names);
   save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials_incorrect', loosestr, '.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    %%
    [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_mid_contrast, parcels_names);
   
    save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials', loosestr, 'mid_contrast.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    
     [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_cor_mid_contrast, parcels_names);
     save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials_correct', loosestr, 'mid_contrast.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    
   [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_inc_mid_contrast, parcels_names);
   save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials_incorrect', loosestr, 'mid_contrast.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    %%
    disp('graph analysis saved')
    if saveplots
        %% Visualization
        % plot graphs
        plot_graph(G_corr, cent_corr_weighted, names_corr, parcels_region_labels);
        suptitle('Correlation');
        set(gcf, 'Position',  [150,150, 1000,500]);
        mysave(gcf, strcat(outputfolder,'correlation_graph',statenames{state_i}), 'fig');
        
        figure;
        for k = 1:length(names_corr)
            subplot(2,4,k);
            bar(cent_corr_weighted.(names_corr{k}));
            set(gca, 'XTickLabel', parcels_names);
            set(gca,'XTickLabelRotation',45);
            set(gcf, 'Position',  [150,150, 1000,500]);
            title(names_corr{k});
        end
        suptitle('Correlation');
        mysave(gcf, strcat(outputfolder,'correlation_centrality',statenames{state_i}), 'fig');
        disp('cntrality plotted and saved')
        
        % plot centrality by node population
        figure;
        subplot(2,1,1);
        bar(sum(indic_corr_weighted,2)/size(indic_corr_weighted, 2));set(gca,'XTickLabel',names_corr);set(gca,'XTickLabelRotation',45);set(gcf, 'Position',  [150,150, 1000,500]);
        title('Correlation');
        
        mysave(gcf, strcat(outputfolder,'centrality_by_node_pop',statenames{state_i}), 'fig');
    end
    
    
end
end
function makeslopeamplitudeplots(animals, isloose,outputfiggolder)
if isloose
    loostr = 'loose';
else
    loostr = '';
end
statenames = {'low_pup_q_zscored', 'high_pup_q_zscored', 'high_pup_l_zscored'};
for state_i=1:length(statenames)
    fields = {'corr','incorr'}; 
    c = cell(length(fields),1);
    x = cell2struct(c,fields);
    for animal_i=1:length(animals)
        animal=char(animals(animal_i));
        res=load(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animals{animal_i},'\',animals{animal_i},'trials_3states', loostr));
        %re = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'), strcat(animaltodays(animal),'imaging_time_traces_global.mat')));
        curr_corr=res.(statenames{state_i}).imaging_time_traces(2,:,res.(statenames{state_i}).trialslabels.blinksummary==1);
        %zscoredbytrial_corr = normByPreActivity(re.imaging_time_traces.t, curr_corr, -3, -1); %normalize activity per day
        curr_incorr=res.(statenames{state_i}).imaging_time_traces(2,:,res.(statenames{state_i}).trialslabels.blinksummary==2);
        
        vectorcorr=mean(squeeze(curr_corr),2);
        vectorincorr=mean(squeeze(curr_incorr),2);
        %Load 10 hz data
        t_10=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\','xz','psych\spt'), strcat(animaltodays('xz'),'imaging_time_traces_global.mat')));
        if strcmp(animal,'xt')||strcmp(animal,'xs')||strcmp(animal,'xu')
            t_33 = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'), strcat(animaltodays(animal),'imaging_time_traces_global.mat')));
            %x_10 = interp1(t_33, x_33, t_10);
            vectorcorr=interp1(t_33.imaging_time_traces.t, vectorcorr, t_10.imaging_time_traces.t).';
            vectorincorr=interp1(t_33.imaging_time_traces.t, vectorincorr, t_10.imaging_time_traces.t).';
        else
        end
        st=findClosestDouble(t_10.imaging_time_traces.t,-0.5);ed=findClosestDouble(t_10.imaging_time_traces.t,2);
        x.corr=cat(2,x.corr,vectorcorr(st:ed));
        x.incorr=cat(2,x.incorr,vectorincorr(st:ed));
        clearvars -except ir animals x t_10 loostr statenames state_i outputfiggolder
    end
t_10=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\','xz','psych\spt'), strcat(animaltodays('xz'),'imaging_time_traces_global.mat')));    
stind=findClosestDouble(t_10.imaging_time_traces.t,-0.5);enind=findClosestDouble(t_10.imaging_time_traces.t,2);
%% plot correct and incorrect
x1 = 0.0; x2 = 0.500;
y1 = -2; y2 = 3;
figure;
set(gcf,'renderer','painters');
fill([x1 x1 x2 x2],[y1 y2 y2 y1],[0.8 0.8 0.8],'LineStyle','none')
hold on
x3 = 0.450; x4 = 0.500;
fill([x3 x3 x4 x4],[y1 y2 y2 y1],[0.3020 0.7490 0.9294],'LineStyle','none')
xlim([-0.5 2])
hold on
shadedErrorBar(transpose(t_10.imaging_time_traces.t(stind:enind)),mean(x.corr,2),(std(x.corr,0,2)./(sqrt(size(x.corr,2)-1))),'lineprops','g');
hold on
shadedErrorBar(transpose(t_10.imaging_time_traces.t(stind:enind)),mean(x.incorr,2),(std(x.incorr,0,2)./(sqrt(size(x.incorr,2)-1))),'lineprops','r');
xlabel('Time [sec]');ylabel('Z-DF/F');title(strcat('Avg Correct vs Incorrect',statenames{state_i}));
hold off
mysave(gcf, fullfile(outputfiggolder,strcat('V1_corr_incorr_',statenames{state_i})), 'all');
clearvars x t_10
end
end
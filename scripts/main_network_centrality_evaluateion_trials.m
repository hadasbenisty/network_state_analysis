function main_network_centrality_evaluateion_trials
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
animals={'xs','xx','xz','xw','xt','xu'};
pre_trial_time_start = -3;
pre_trial_time_end = -.1;
isloose = true;
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
for ai = 1:length(animals)
%     eval_weights_and_cent(isloose, animals{ai}, statenames, pre_trial_time_start, pre_trial_time_end);
end
ismidcontrast=false;
outputfiggolder = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
%% gal

plotSummaryCentrality_gal(isloose, animals, outputfiggolder, statenames);

%% allen
[trials_states_notweighted, trials_states_weighted] = plot_centrality_res(ismidcontrast, isloose, animals, outputfiggolder, statenames, 'trials');
[correct_states_notweighted, correct_states_weighted] = plot_centrality_res(ismidcontrast, isloose, animals, outputfiggolder, statenames, 'trials_correct');
[incorrect_states_notweighted, incorrect_states_weighted] = plot_centrality_res(ismidcontrast, isloose, animals, outputfiggolder, statenames, 'trials_incorrect');
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
save(fullfile(outputfiggolder, ['centrality_stats_pretrial' loosestr midcontraststr '.mat']), 'trials_states_notweighted',...
    'trials_states_weighted', 'correct_states_notweighted', 'correct_states_weighted',...
    'incorrect_states_notweighted', 'incorrect_states_weighted');
plotSummaryCentrality(ismidcontrast, isloose, outputfiggolder, statenames);
makeslopeamplitudeplots(animals, isloose,outputfiggolder)

end
function plotSummaryCentrality_gal(isloose, animals, outputfiggolder, statenames)

[correct_states_notweighted] = plot_centrality_res_gal(isloose, animals, outputfiggolder, statenames, 'trials_correct');
[incorrect_states_notweighted] = plot_centrality_res_gal(isloose, animals, outputfiggolder, statenames, 'trials_incorrect');

cent_features = fieldnames(correct_states_notweighted.high_pup_l);
    [parcels_names, ~, finalindex] = get_allen_meta_parcels;
parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
for state_i=1:length(statenames)
for ni = 1:length(cent_features)
    figure;
    
    %difference maps, not weighted, for each centrality measure
     Ac=mean(correct_states_notweighted.(statenames{state_i}).(cent_features{ni}),3);
     Ain=mean(incorrect_states_notweighted.(statenames{state_i}).(cent_features{ni}),3);
     A=Ac-Ain;
    imagesc(A);
set(gcf,'renderer','painters');
myColorMap = colormap(redblue);
%myColorMap(1,:) = 1;
colormap(myColorMap);
h=colorbar;
upperlim=max(A(:)); 
    lowerlim=min(A(:));
caxis([lowerlim upperlim]);
%colormap(fireice);h=colorbar;
ylabel(h, 'Difference in Node Centrality');
title([cent_features{ni} ' ' statenames{state_i}]);axis off
hold on
plot_parcellation_boundaries(parcelsallen.parcells_new.indicators(:,:,finalindex));

end
   
  
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
    for k=1:length(statenames)
        M(:, 1, k) = nanmean(correct_states_notweighted.(statenames{k}).(centstr),2);
        M(:, 2, k) = nanmean(incorrect_states_notweighted.(statenames{k}).(centstr),2);
        S(:, 1, k) = nanstd(correct_states_notweighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
        S(:, 2, k) = nanstd(incorrect_states_notweighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
    end
    plot_correct_incorrect_per_state_per_parcels(M, S, parcels_names, statenames)
    mysave(gcf, fullfile(outputfiggolder, 'not_weighted', [centnames{l} '_centrality_stats_pretrial' loosestr midcontraststr]));
end


% for l=1:length(centnames)
%     centstr = centnames{l};
%     figure;
% for k=1:length(statenames)
% M1 = nanmean(correct_states_notweighted.(statenames{k}).(centstr),2);
% S1 = nanstd(correct_states_notweighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
% M2 = nanmean(incorrect_states_notweighted.(statenames{k}).(centstr),2);
% S2 = nanstd(incorrect_states_notweighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
%
% subplot(3, 1, k);
% barwitherr([S1 S2], [M1 M2]);
% set(gca,'XTick', 1:length(parcels_names))
% set(gca,'XTickLabel', parcels_names)
% strttl = statenames{k};
% strttl(strttl=='_') = ' ';
% title(strttl);
% axis tight;
% end
% legend('Correct','Incorrect');
% suptitle(['Correct/Incorrect ' centstr]);
% set(gcf, 'Position', [1          41        1920         963]);
%
%
% mysave(gcf, fullfile(outputfiggolder, 'not_weighted', [centnames{l} '_centrality_stats_pretrial' loosestr midcontraststr]));
% end
% %%  weighted
% centnames = fieldnames(correct_states_weighted.low_pup_q);
% for l=1:length(centnames)
%     centstr = centnames{l};
%     figure;
% for k=1:length(statenames)
% M1 = nanmean(correct_states_weighted.(statenames{k}).(centstr),2);
% S1 = nanstd(correct_states_weighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
% M2 = nanmean(incorrect_states_weighted.(statenames{k}).(centstr),2);
% S2 = nanstd(incorrect_states_weighted.(statenames{k}).(centstr),[],2)/sqrt(N-1);
% subplot(3, 1, k);
% barwitherr([S1 S2], [M1 M2]);
% set(gca,'XTick', 1:length(parcels_names))
% set(gca,'XTickLabel', parcels_names)
% strttl = statenames{k};
% strttl(strttl=='_') = ' ';
% title(strttl);
% axis tight;
% end
% legend('Correct','Incorrect');
% suptitle(['Correct/Incorrect ' centstr]);
% set(gcf, 'Position', [1          41        1920         963]);
% mysave(gcf, fullfile(outputfiggolder, 'weighted', [centnames{l} '_centrality_stats_pretrial' loosestr  midcontraststr]));
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
%%
function [spon_states_notweighted, spon_states_weighted] = plot_centrality_res_gal(isloose, animals, outputfiggolder, statenames, suffix_files)
if isloose
    loosestr = 'loose';
else
    loosestr = '';
end

[parcels_names] = get_allen_meta_parcels;
cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector', 'participation'};
for state_i = 1:length(statenames)
    for cent_i = 1:length(cent_features)
        spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = zeros(256, 256, length(animals));
        spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = zeros(256, 256, length(animals));
    end
end
for i=1:length(animals)
    animal=animals{i};
    [parcels_names, ~, finalindex] = get_allen_meta_parcels;
    [roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal] = get_gal_parcels_lables(animal);

[parcels_names_gal, finalindex_gal, maskByGal, regionLabel_gal] = get_gal_meta_parcels_by_allen(parcels_names, finalindex,...
    roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal);

    outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020',animal);
    for state_i = 1:length(statenames)
        
        
        load(fullfile(outputfolder,['network_analysis_corr',statenames{state_i} ,suffix_files, loosestr, 'gal.mat']), ...
            'cent_corr_weighted','cent_corr_notweighted');
        cent_features = fieldnames(cent_corr_weighted);
        for cent_i = 1:length(cent_features)
            P=zeros(256);
            for parcel_i = 1:length(parcels_names_gal)
                 P(maskByGal==parcel_i) = cent_corr_notweighted.(cent_features{cent_i})(parcel_i);
            end
            spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i})(:,:,i) = P;               
            
        end
        
    end
end
mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
legstr = {'Low Q', 'High Q', 'Loc'};
 parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
%   
% for ni = 1:length(cent_features)
%     figure;
%     
%     %difference maps, not weighted, for each centrality measure
%      A1=mean(spon_states_notweighted.(statenames{1}).(cent_features{ni}),3);
%      A3=mean(spon_states_notweighted.(statenames{3}).(cent_features{ni}),3);
%      A=A3-A1;
%     imagesc(A);
% set(gcf,'renderer','painters');
% myColorMap = colormap(redblue);
% %myColorMap(1,:) = 1;
% colormap(myColorMap);
% h=colorbar;
% upperlim=max(A(:)); 
%     lowerlim=min(A(:));
% caxis([lowerlim upperlim]);
% %colormap(fireice);h=colorbar;
% ylabel(h, 'Difference in Node Centrality');title(cent_features{ni});axis off
% hold on
% plot_parcellation_boundaries(parcelsallen.parcells_new.indicators(:,:,finalindex));
% 
%  
%    
%   
% end

end

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

function eval_weights_and_cent(isloose, animal, statenames, pre_trial_time_start, pre_trial_time_end)
if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
mkNewDir(outputfolder);
[roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal] = get_gal_parcels_lables(animal);
load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states' loosestr '.mat'],...
    'low_pup_q','high_pup_q','high_pup_l','days_to_process', 't'); %#ok<NASGU>
[parcels_names, ~, finalindex] = get_allen_meta_parcels;
[parcels_names_gal, finalindex_gal, maskByGal, regionLabel_gal] = get_gal_meta_parcels_by_allen(parcels_names, finalindex,...
    roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal);

disp(animal)
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    
    data_3D = eval(statenames{state_i});
    imaging_time_traces_all = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary<3);
    imaging_time_traces_cor = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary==1);
    imaging_time_traces_inc = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary==2);
    
    imaging_time_traces_all_gal = data_3D.Gal(:, :, data_3D.trialslabels.blinksummary<3);
    imaging_time_traces_cor_gal = data_3D.Gal(:, :, data_3D.trialslabels.blinksummary==1);
    imaging_time_traces_inc_gal = data_3D.Gal(:, :, data_3D.trialslabels.blinksummary==2);
    
    
    data=[];data_inco=[];data_corr=[];
    data_gal=[];data_inco_gal=[];data_corr_gal=[];
    for T=1:size(imaging_time_traces_all,3)
        data = cat(2, data,  imaging_time_traces_all(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_cor,3)
        data_corr = cat(2, data_corr,  imaging_time_traces_cor(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_inc,3)
        data_inco = cat(2, data_inco,  imaging_time_traces_inc(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_all_gal,3)
        data_gal = cat(2, data_gal,  imaging_time_traces_all_gal(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_cor_gal,3)
        data_corr_gal = cat(2, data_corr_gal,  imaging_time_traces_cor_gal(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    for T=1:size(imaging_time_traces_inc_gal,3)
        data_inco_gal = cat(2, data_inco_gal,  imaging_time_traces_inc_gal(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
    end
    
    data = data(finalindex, :);
    data_corr = data_corr(finalindex, :);
    data_inco = data_inco(finalindex, :);
    
    data_gal = data_gal(finalindex_gal, :);
    data_corr_gal = data_corr_gal(finalindex_gal, :);
    data_inco_gal = data_inco_gal(finalindex_gal, :);
    
    data = data(:, all(~isnan(data)));
    data_corr = data_corr(:, all(~isnan(data_corr)));
    data_inco = data_inco(:, all(~isnan(data_inco)));
    
    data_gal = data_gal(:, all(~isnan(data_gal)));
    data_corr_gal = data_corr_gal(:, all(~isnan(data_corr_gal)));
    data_inco_gal = data_inco_gal(:, all(~isnan(data_inco_gal)));
    if any(isnan(data(:)))
        disp('nans in dataset')
        continue;
    end
    %% measure weights
    W_corr = measure_weights_partial(data, 'corr');
    W_corr_cor = measure_weights_partial(data_corr, 'corr');
    W_corr_inc = measure_weights_partial(data_inco, 'corr');
    
    W_corr_gal = measure_weights_partial(data_gal, 'corr');
    W_corr_cor_gal = measure_weights_partial(data_corr_gal, 'corr');
    W_corr_inc_gal = measure_weights_partial(data_inco_gal, 'corr');
    %
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
    %% gal
    [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_gal, parcels_names);
    W_corr=W_corr_gal;
    save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials', loosestr, 'gal.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    
    [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_cor_gal, parcels_names);
    W_corr=W_corr_cor_gal;
    save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials_correct', loosestr, 'gal.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    
    [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_inc_gal, parcels_names);
    W_corr=W_corr_inc_gal;
    save(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'trials_incorrect', loosestr, 'gal.mat'),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    %%
    disp('graph analysis saved')
    
    
    
end
end

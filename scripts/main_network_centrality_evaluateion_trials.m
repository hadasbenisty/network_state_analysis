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
for ai = 1%:length(animals)
     %eval_weights_and_cent(isloose, animals{ai}, statenames, pre_trial_time_start, pre_trial_time_end);
%      eval_weights_and_cent_perm(isloose, animals{ai}, statenames,pre_trial_time_start, pre_trial_time_end)
end
ismidcontrast=false;
outputfiggolder = 'X:\Lav\ProcessingDirectory\parcor_undirected\';

%% gal

%plotSummaryCentrality_gal(isloose, animals, outputfiggolder, statenames);

%% allen
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
% plotSummaryCentrality(ismidcontrast, isloose, outputfiggolder, statenames);
 plotSummaryCentrality_gal(isloose, animals, outputfiggolder, statenames)
[~,spatialindex]=getspatialindex;
makeslopeamplitudeplots(animals, isloose,outputfiggolder,spatialindex(1),'V1')
makeslopeamplitudeplots(animals, isloose,outputfiggolder,spatialindex(2),'S1b')
makeslopeamplitudeplots(animals, isloose,outputfiggolder,spatialindex(3),'M2')

%plotCRF(isloose,animals)
end
function plotSummaryCentrality_gal(isloose, animals, outputfiggolder, statenames)
if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
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
mysave(gcf, fullfile(outputfiggolder, 'not_weighted', 'trials',[statenames{state_i} cent_features{ni} '_centrality_stats_pretrial_gal' loosestr '_heatmap']));

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
function eval_weights_and_cent_perm(isloose, animal, statenames,pre_trial_time_start, pre_trial_time_end)
if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
mkNewDir(outputfolder);
load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states' loosestr '.mat'],...
    'low_pup_q','high_pup_q','high_pup_l','days_to_process', 't'); %#ok<NASGU>
[parcels_names, ~, finalindex] = get_allen_meta_parcels;

disp(animal)
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    
    data_3D = eval(statenames{state_i});

    imaging_time_traces_cor = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary==1);
    imaging_time_traces_inc = data_3D.imaging_time_traces(:, :, data_3D.trialslabels.blinksummary==2);
    perm_data_corr=[];perm_data_incorr=[];
    for perm_i=1:100
        [corralldata,incorralldata]=permute_corrincorr_lan(imaging_time_traces_cor,imaging_time_traces_inc);
        for T=1:size(corralldata,3)
            perm_data_corr = cat(2, perm_data_corr,  corralldata(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
        end
        for T=1:size(incorralldata,3)
            perm_data_incorr = cat(2, perm_data_incorr,  incorralldata(:,t>=pre_trial_time_start & t<pre_trial_time_end,T));
        end
    save(strcat(outputfolder,'\perm\',num2str(perm_i),statenames{state_i},'corrincorr_perm'),'perm_data_corr','perm_data_incorr')
    data=load(strcat(outputfolder,'\perm\',num2str(perm_i),statenames{state_i},'corrincorr_perm'));

    %change loading and saving
    data_corr = data.perm_data_corr(:, all(~isnan(data.perm_data_corr)));
    data_inco = data.perm_data_incorr(:, all(~isnan(data.perm_data_incorr)));
   
    if any(isnan(data_corr(:)))|any(isnan(data_inco(:)))
        disp('nans in dataset')
        continue;
    end
    %% measure weights
    W_corr_cor = measure_weights_partial(data_corr, 'corr');
    W_corr_inc = measure_weights_partial(data_inco, 'corr');
    
    disp('corweights done')
    
    [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_cor, parcels_names);
    save(strcat(outputfolder,num2str(perm_i),'network_analysis_corr_perm',statenames{state_i} ,'trials_correct', loosestr, '.mat'),'W_corr_cor',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    
    [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr_inc, parcels_names);
    save(strcat(outputfolder,num2str(perm_i),'network_analysis_corr_perm',statenames{state_i} ,'trials_incorrect', loosestr, '.mat'),'W_corr_inc',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    disp('graph analysis saved')
    end   
end
end
    
    
    
function makeslopeamplitudeplots(animals, isloose,outputfiggolder,spatialindex,parcelname)

if isloose

    loostr = 'loose';

else

    loostr = '';
end
cc=20;
statenames = {'low_pup_q_zscored', 'high_pup_q_zscored', 'high_pup_l_zscored'};
for state_i=1:length(statenames)
    fields = {'corr','incorr'}; 
    c = cell(length(fields),1);
    x = cell2struct(c,fields);
    for animal_i=1:length(animals)
        animal=char(animals(animal_i));
        res=load(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animals{animal_i},'\',animals{animal_i},'trials_3states', loostr));
        %re = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'), strcat(animaltodays(animal),'imaging_time_traces_global.mat')));
        indxcorr=res.(statenames{state_i}).trialslabels.blinksummary==1&res.(statenames{state_i}).trialslabels.contrastLabels==cc;
        indxincorr=res.(statenames{state_i}).trialslabels.blinksummary==2&res.(statenames{state_i}).trialslabels.contrastLabels==cc;

        curr_corr=res.(statenames{state_i}).imaging_time_traces(spatialindex,:,indxcorr);
        %zscoredbytrial_corr = normByPreActivity(re.imaging_time_traces.t, curr_corr, -3, -1); %normalize activity per day
        curr_incorr=res.(statenames{state_i}).imaging_time_traces(spatialindex,:,indxincorr);

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
        clearvars -except cc parcelname ir animals x t_10 loostr statenames state_i outputfiggolder spatialindex
    end
t_10=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\','xz','psych\spt'), strcat(animaltodays('xz'),'imaging_time_traces_global.mat')));    
stind=findClosestDouble(t_10.imaging_time_traces.t,-0.5);enind=findClosestDouble(t_10.imaging_time_traces.t,2);
%% plot correct and incorrect
x1 = 0.0; x2 = 0.500;
y1 = -2; y2 = 12;
figure;
set(gcf,'renderer','painters');
fill([x1 x1 x2 x2],[y1 y2 y2 y1],[0.8 0.8 0.8],'LineStyle','none')
hold on
x3 = 0.450; x4 = 0.500;
fill([x3 x3 x4 x4],[y1 y2 y2 y1],[0.3020 0.7490 0.9294],'LineStyle','none')
xlim([-0.5 2]);ylim([-2 5])
hold on
shadedErrorBar(transpose(t_10.imaging_time_traces.t(stind:enind)),mean(x.corr,2),(std(x.corr,0,2)./(sqrt(size(x.corr,2)-1))),'lineprops','g');
hold on
shadedErrorBar(transpose(t_10.imaging_time_traces.t(stind:enind)),mean(x.incorr,2),(std(x.incorr,0,2)./(sqrt(size(x.incorr,2)-1))),'lineprops','r');
xlabel('Time [sec]');ylabel('Z-DF/F');title(strcat('Avg Correct vs Incorrect',statenames{state_i}));
hold off
mysave(gcf, fullfile(outputfiggolder,strcat(parcelname,num2str(cc),'contrast_corr_incorr_',statenames{state_i})), 'all');
clearvars x t_10
end
end

function plotCRF(isloose,animals)
if isloose
    loostr = 'loose';
else
    loostr = '';
end
contrasts=[0 2 5 10 20 40 100];
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
acs=NaN(length(animals),length(statenames),length(contrasts));
%option for subtracting 300ms prior to stim or no normalization. currently i am saving
%zscored to mean and std in 300 ms prior to stim
acs_sub=NaN(length(animals),length(statenames),length(contrasts));
acs_z=NaN(length(animals),length(statenames),length(contrasts));
for animal_i=1:length(animals)
    res=load(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animals{animal_i},'\',animals{animal_i},'trials_3states', loostr));
    for state_i=1:length(statenames)
        %concatenate mean response per animal in 300 ms
        for contrast_i=1:length(contrasts)
            indx=res.(statenames{state_i}).trialslabels.contrastLabels==contrasts(contrast_i)&res.(statenames{state_i}).trialslabels.blinksummary<3;
            v1_act=squeeze(res.(statenames{state_i}).imaging_time_traces(2,:,indx));
            %row wise zscoring
            st_z=findClosestDouble(res.t,-3);ed_z=findClosestDouble(res.t,0);            
            v1_act_norm1 = bsxfun(@minus, v1_act, nanmean(v1_act(st_z:ed_z,:)));
            v1_act_norm2= bsxfun(@rdivide, v1_act_norm1, nanstd(v1_act(st_z:ed_z,:)));
            st=findClosestDouble(res.t,0);
            ed=findClosestDouble(res.t,0.3);
            acs(animal_i,state_i,contrast_i)=nanmean(nanmean(v1_act(st:ed,:),1));
            acs_sub(animal_i,state_i,contrast_i)=nanmean(nanmean(v1_act_norm1(st:ed,:),1));
            acs_z(animal_i,state_i,contrast_i)=nanmean(nanmean(v1_act_norm2(st:ed,:),1));
            clearvars indx v1_act st ed
        end
    end
end
CondColors=[0,0,0;0.9290 0.6940 0.1250;1,0,0];
figure;
%hAx=axes;
%hAx.XScale='log';
xlim([0 100]);
hold on
for state_i=1:length(statenames)
    M=squeeze(acs_z(:,state_i,:));
    %errorbar(contrasts,nanmean(M,1),nanstd(M,1)./(sqrt(length(animals)-1)));
    %[modelF,coefsF]=HyFit(contrasts.',nanmean(M,1));
    %fittedpredictedcurveF=nanmean(predint(modelF,0:0.01:100),2); 
    %semilogx(0:0.01:100,fittedpredictedcurveF,'color',CondColors(state_i,:)); %overlay early and late days psyc curve on semi log plo
    errorbar(contrasts,nanmean(M,1),nanstd(M,1)./(sqrt(length(animals)-1)),'color',CondColors(state_i,:),'LineWidth',1.5);
    hold on
end
New_XTickLabel = get(gca,'xtick');title('CRF per state');ylabel('Z-dF/F')
set(gca,'XTickLabel',New_XTickLabel);
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\CRF_Zscored_perstate', 'all');
end

function [corralldata,incorralldata]=permute_corrincorr_lan(imaging_time_traces_cor,imaging_time_traces_inc)
%concatenate conditions
concatenated_dff=cat(3,imaging_time_traces_cor,imaging_time_traces_inc);
%size of each condition and size of the entire dataset to permute
samplesize=size(imaging_time_traces_cor,3);
permutation = randperm(size(concatenated_dff,3));
%save permutations
corralldata=concatenated_dff(:,:,permutation(1:samplesize));
incorralldata=concatenated_dff(:,:,permutation(samplesize+1:length(permutation)));
%save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',num2str(i),animal,'perm_running_time_traces'),'runningalldata','runningalldata_t','notrunningalldata_t','notrunningalldata');
end


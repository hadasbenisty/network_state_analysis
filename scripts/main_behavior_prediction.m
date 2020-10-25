function main_behavior_prediction
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../svm'));
animals={'xs','xx','xz','xw','xt','xu'};
slope_trial_time_start = 0.1;
slope_trial_time_end = 0.33;
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
for ai = 1:length(animals)
    behavior_prediction(animals{ai}, statenames, slope_trial_time_start, slope_trial_time_end);
end
outputfiggolder = 'X:\Lav\ProcessingDirectory\parcor_undirected\';


plot_prediction(animals, statenames, outputfiggolder);
end
function plot_prediction(animals, statenames, outputfiggolder)


for ai = 1:length(animals)
    
    outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animals{ai},'\');
    for state_i = 1:length(statenames)
        if ~exist(strcat(outputfolder,'behavior_prediction',statenames{state_i} ,'.mat'), 'file')
            acc_all_parcels(state_i, ai) = nan;
            acc_per_parcel(:, state_i, ai) = nan(23,1);
        else
            load(strcat(outputfolder,'behavior_prediction',statenames{state_i} ,'.mat'),'slopeData',...
                'accuracy_vec','cvinds','fa_vec',...
                'md_vec', 'accuracy_mat', 'cvinds_mat', 'fa_mat', 'md_mat');
            acc_all_parcels(state_i, ai) = mean(accuracy_vec);
            acc_per_parcel(:, state_i, ai) = mean(accuracy_mat);
        end
    end
end
n=length(animals);
M = nanmean(acc_all_parcels,2);
S = nanstd(acc_all_parcels,[],2)/sqrt(n-1);
figure;barwitherr(S,M);
set(gca,'XTickLabel', statenames);
ylim([0.5 0.9]);
mysave(gcf,fullfile(outputfiggolder, 'behavior_prediction_by_state_all_parcels'));
parcels_names = get_allen_meta_parcels;
figure;
M = nanmean(acc_per_parcel,3);
S = nanstd(acc_per_parcel,[],3)/sqrt(n-1);
barwitherr(S,M);
legend(statenames)
ylim([0.5 0.9]);
set(gca,'XTick', 1:length(parcels_names));
set(gca,'XTickLabel', parcels_names);
set(gcf,'Position',[1          41        1920         963])
mysave(gcf,fullfile(outputfiggolder, 'behavior_prediction_by_state_per_parcel'));
end
function behavior_prediction(animal, statenames, slope_trial_time_start, slope_trial_time_end)
outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
mkNewDir(outputfolder);

[~, ~, finalindex] = get_allen_meta_parcels;
load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states.mat'],...
    'low_pup_q','high_pup_q','high_pup_l','days_to_process', 't'); %#ok<NASGU>
disp(animal);
foldsNum=5;
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    
    data_3D = eval(statenames{state_i});
    inds = data_3D.trialslabels.blinksummary<3;
    imaging_time_traces = data_3D.imaging_time_traces(:, :, inds);
    labels = data_3D.trialslabels.blinksummary(inds);
    imaging_time_traces = imaging_time_traces(finalindex, :, :);
    slopeData=[];accuracy_mat=[];cvinds_mat=[];fa_mat=[];md_mat=[];
    for T=1:size(imaging_time_traces,3)
        for pari = 1:size(imaging_time_traces,1)
        p = polyfit(t(t>=slope_trial_time_start & t<slope_trial_time_end),imaging_time_traces(pari,t>=slope_trial_time_start & t<slope_trial_time_end,T),1);
        slopeData(pari, T) = p(1);
        end
    end
    [accuracy_vec, cvinds, fa_vec, md_vec]  = svmClassify(slopeData.', labels, foldsNum, 0, 1, 'downsample', 0, 0, 0, 0);   
   if isempty(accuracy_vec)
       continue;
   end
    for par_i = 1:size(imaging_time_traces,1)
    [accuracy_mat(:, par_i), cvinds_mat(:, par_i), fa_mat(:,:, par_i), md_mat(:, :,par_i)]  = svmClassify(slopeData(par_i,:).', labels, foldsNum, 0, 1, 'downsample', 0, 0, 0, 0);   
   end
    
    save(strcat(outputfolder,'behavior_prediction',statenames{state_i} ,'.mat'),'slopeData',...
        'accuracy_vec','cvinds','fa_vec',...
        'md_vec', 'accuracy_mat', 'cvinds_mat', 'fa_mat', 'md_mat');
    
    
    
end
end






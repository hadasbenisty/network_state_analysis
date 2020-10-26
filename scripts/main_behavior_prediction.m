function main_behavior_prediction
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../svm'));
animals={'xs','xx','xz','xw','xt','xu'};
slope_trial_time_start = 0.1;
slope_trial_time_end = 0.33;
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
isloose = true;
for ai = 1:length(animals)
    behavior_prediction(isloose, animals{ai}, statenames, slope_trial_time_start, slope_trial_time_end);
end
outputfiggolder = 'X:\Lav\ProcessingDirectory\parcor_undirected\';

contrast_levels = [0 2 5 10 20 40 100];

% plot_prediction(isloose, animals, statenames, outputfiggolder);
plot_psych_curve_per_state(isloose, animals, statenames, contrast_levels, outputfiggolder)
end

function plot_psych_curve_per_state(isloose, animals, statenames, contrast_levels, outputfiggolder)

if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
suc_rate=nan(length(statenames), length(animals));
psych_curv=nan(length(contrast_levels), length(statenames), length(animals));
for ai=1:length(animals)

load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animals{ai},'\',animals{ai},'trials_3states', loosestr, '.mat'],...
    'low_pup_q','high_pup_q','high_pup_l'); %#ok<NASGU>

for state_i = 1:length(statenames)
    disp(statenames{state_i})
    
    data_3D = eval(statenames{state_i});
    inds = data_3D.trialslabels.blinksummary<3;    
    labels = data_3D.trialslabels.blinksummary(inds);
    contrast = data_3D.trialslabels.contrastLabels(inds);    
    suc_rate(state_i, ai) = sum(labels==1)/length(labels);
    for ci=1:length(contrast_levels)
        psych_curv(ci, state_i, ai) =  sum(labels==1 & contrast == contrast_levels(ci))/sum(contrast == contrast_levels(ci));        
    end
end
end
n = length(animals); 
M = nanmean(psych_curv, 3)*100;
S = nanstd(psych_curv, [], 3)/sqrt(n-1)*100;

figure;b=errorbar(contrast_levels, M(:,1), S(:,1), 'b');
hold all;
errorbar(contrast_levels, M(:,2), S(:,2), 'Color',[0.9290 0.6940 0.1250]);
errorbar(contrast_levels, M(:,3), S(:,3),'r');
xlabel('% Contrast');
ylabel('% Success Rate');
legend(statenames);
mysave(gcf, fullfile(outputfiggolder,['psych_curves_by_state' loosestr]));
  
n = length(animals); 
M = nanmean(suc_rate, 2)*100;
S = nanstd(suc_rate, [], 2)/sqrt(n-1)*100;
figure;barwitherr(S, M)
set(gca,'XTickLabel',statenames);
ylabel('% Success Rate');
mysave(gcf, fullfile(outputfiggolder,['success_rate_by_state' loosestr]));

    

end


function plot_prediction(isloose, animals, statenames, outputfiggolder)

if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
for ai = 1:length(animals)
    
    outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animals{ai},'\');
    for state_i = 1:length(statenames)
        if ~exist(strcat(outputfolder,'behavior_prediction',statenames{state_i}, loosestr ,'.mat'), 'file')
            acc_all_parcels(state_i, ai) = nan;
            acc_per_parcel(:, state_i, ai) = nan(23,1);
        else
            load(strcat(outputfolder,'behavior_prediction',statenames{state_i}, loosestr ,'.mat'),'slopeData',...
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
mysave(gcf,fullfile(outputfiggolder, ['behavior_prediction_by_state_all_parcels' loosestr]));
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
mysave(gcf,fullfile(outputfiggolder, ['behavior_prediction_by_state_per_parcel' loosestr]));
end
function behavior_prediction(isloose, animal, statenames, slope_trial_time_start, slope_trial_time_end)

if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
mkNewDir(outputfolder);

[~, ~, finalindex] = get_allen_meta_parcels;
load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states' loosestr '.mat'],...
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
    
    save(strcat(outputfolder,'behavior_prediction',statenames{state_i}, loosestr,'.mat'),'slopeData',...
        'accuracy_vec','cvinds','fa_vec',...
        'md_vec', 'accuracy_mat', 'cvinds_mat', 'fa_mat', 'md_mat');
    
    
    
end
end






% Figure 3
function main_behavior_prediction_2p_lan
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../svm'));
T = readtable('X:\Lan\FMB208 conditioning imaging data\animal_log_imaging_lt.xlsx');
animals = T.AnimalID;

slope_trial_time_start = 0.1;
slope_trial_time_end = 0.33;
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
daysmax=30;
for ai = 1:length(animals)
    behavior_prediction(animals{ai}, statenames, slope_trial_time_start, slope_trial_time_end,daysmax);
end
end
function behavior_prediction(animal, statenames, slope_trial_time_start, slope_trial_time_end,daysmax)


outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\behavior_prediction';
mkNewDir(outputfolder);
metafile = ['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_meta_data.mat'];
imaginggile = ['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_imaging_time_traces_global_dfff.mat'];
if ~exist(metafile, 'file') || ~exist(imaginggile, 'file')
    return;
end
load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_meta_data.mat'], ...
    'behavior_labels', 'arousal_labels', 'days_labels');

load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_imaging_time_traces_global_dfff.mat'], ...
    'imaging_time_tracesN', 't' );
       
disp(animal);
foldsNum=5;
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    if exist(fullfile(outputfolder,[animal '_' statenames{state_i} '.mat']),'file')
        continue;
    end
   
    inds = behavior_labels<3 & arousal_labels == state_i&days_labels<=daysmax;
    curr_imaging_time_traces = imaging_time_tracesN(:, :, inds);
    labels = behavior_labels(inds);
    
    slopeData=[];accuracy_mat=[];
    for T=1:size(curr_imaging_time_traces,3)        
        slopeData(:,T) = getSlope(curr_imaging_time_traces(:, t>=slope_trial_time_start & t<slope_trial_time_end, T));
    end
    accuracy_vec  = svmClassify(slopeData.', labels, foldsNum, 0, 1, 'downsample', 0, 0, 0, 0);
    if isempty(accuracy_vec)
        continue;
    end
      for par_i = 1:size(slopeData,1)
        accuracy_mat(:, par_i)  = svmClassify(slopeData(par_i,:).', labels, foldsNum, 0, 1, 'downsample', 0, 0, 0, 0);
    end
  
    
    save(fullfile(outputfolder,[animal '_' statenames{state_i} '.mat']),'slopeData',...
        'accuracy_vec',...
         'accuracy_mat');
end

end

























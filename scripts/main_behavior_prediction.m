function main_behavior_prediction
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../svm'));
animals={'xw','xx','xt','xu' 'xs','xz',};
slope_trial_time_start = 0.1;
slope_trial_time_end = 0.33;
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};

for ai = 1:length(animals)
    behavior_prediction(animals{ai}, statenames, slope_trial_time_start, slope_trial_time_end);
end
end
function behavior_prediction(animal, statenames, slope_trial_time_start, slope_trial_time_end)


outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\behavior_prediction';
mkNewDir(outputfolder);
signalsnames = {'Allen','LSSC','Grid4'};
load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_meta_data.mat'], ...
    'behavior_labels', 'arousal_labels');

for sig_i = 1:length(signalsnames)
load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '_trial_imaging_time_traces_global_' signalsnames{sig_i} '.mat'], ...
    'imaging_time_traces', 't' , 'roi_location', 'region_labels', 'roi_names');
  
disp(animal);
foldsNum=5;
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    
   
    inds = behavior_labels<3 & arousal_labels == state_i;
    curr_imaging_time_traces = imaging_time_traces(:, :, inds);
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
  
    
    save(fullfile(outputfolder,[animal '_' statenames{state_i}, '_' signalsnames{sig_i} '.mat']),'slopeData',...
        'accuracy_vec',...
         'accuracy_mat', 'roi_location', 'region_labels', 'roi_names');
end
end
end
function slopeData = getSlope(X)

for pari = 1:size(X,1)
    p = polyfit(1:size(X,2),X(pari,:),1);
    slopeData(pari) = p(1); %#ok<AGROW>
end
end
























% 
% 
% function behavior_prediction_grid(isloose, animal, statenames, slope_trial_time_start, slope_trial_time_end, Gv)
% 
% if isloose
%     loosestr = 'loose';
% else
%     loosestr = '';
% end
% outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
% mkNewDir(outputfolder);
% [parcels_names, ~, finalindex] = get_allen_meta_parcels;
% 
% 
% for g=1:length(Gv)
%    
%     load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states_grid' num2str(Gv(g)) loosestr '.mat'],...
%         'low_pup_q','high_pup_q','high_pup_l','days_to_process', 't', 'parinds');
%      [parcels_names_grid, parcels_region_labels_grid, final_index_grid, region_lut,...
%         grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(parinds, Gv(g));
%     disp(animal);
%     foldsNum=5;
%     for state_i = 1:length(statenames)
%         disp(statenames{state_i})
%         
%         data_3D = eval(statenames{state_i});
%         inds = data_3D.trialslabels.blinksummary<3;
%         imaging_time_traces = data_3D.grid(:, :, inds);
%         labels = data_3D.trialslabels.blinksummary(inds);
%         imaging_time_traces = imaging_time_traces(final_index_grid, :, :);
%         slopeData=[];accuracy_mat=[];cvinds_mat=[];fa_mat=[];md_mat=[];
%         
%         for T=1:size(imaging_time_traces,3)
%             for pari = 1:size(imaging_time_traces,1)
%                 p = polyfit(t(t>=slope_trial_time_start & t<slope_trial_time_end),imaging_time_traces(pari,t>=slope_trial_time_start & t<slope_trial_time_end,T),1);
%                 slopeData(pari, T) = p(1);
%             end           
%         end
%         [accuracy_vec, cvinds, fa_vec, md_vec]  = svmClassify(slopeData.', labels, foldsNum, 0, 1, 'downsample', 0, 0, 0, 0);
%         if isempty(accuracy_vec)
%             continue;
%         end
%         
%         for par_i = 1:size(imaging_time_traces,1)
%             [accuracy_mat(:, par_i), cvinds_mat(:, par_i), fa_mat(:,:, par_i), md_mat(:, :,par_i)]  = svmClassify(slopeData(par_i,:).', labels, foldsNum, 0, 1, 'downsample', 0, 0, 0, 0);
%         end
%       
%         save(strcat(outputfolder,'behavior_prediction',statenames{state_i}, '_grid', num2str(Gv(g)), loosestr,'.mat'),'slopeData',...
%             'accuracy_vec','cvinds','fa_vec',...
%             'md_vec', 'accuracy_mat', 'cvinds_mat', 'fa_mat', 'md_mat');
%         c=unique(grid_map_final_index);
%         A=zeros(size(grid_map_final_index));
%         v=mean(accuracy_mat);
%         for i=2:length(c)
%             A(grid_map_final_index==c(i)) = v(i-1);
%             
%         end
%         figure;imagesc(A,[0.5 .8]);colormap(redblue)
% title(statenames{state_i})
% 
%         clear accuracy_mat;clear md_vec;clear fa_mat;
%         clear cvinds;clear fa_vec;clear cvinds_mat;clear md_mat;
%     end
% end
% end
% 
% 

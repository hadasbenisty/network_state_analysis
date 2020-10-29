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
%     behavior_prediction(isloose, animals{ai}, statenames, slope_trial_time_start, slope_trial_time_end);
end
outputfiggolder = 'X:\Lav\ProcessingDirectory\parcor_undirected\';

contrast_levels = [0 2 5 10 20 40 100];

plot_prediction(isloose, animals, statenames, outputfiggolder);
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
%% c50 and rmax per animal for barplots
c50s=NaN(length(animals),length(statenames));
rmax=NaN(length(animals),length(statenames));
baseline=NaN(length(animals),length(statenames));
for animal_i=1:length(animals)
    for state_i=1:length(statenames)
        [~,coefs]=HyFit([0 2 5 10 20 40 100].',psych_curv(:,state_i,animal_i));
        c50s(animal_i,state_i)=coefs(3);
        baseline(animal_i,state_i)=coefs(4);
        rmax(animal_i,state_i)=coefs(1);
        clearvars coefs
    end
end

makepsychbarplot(c50s,animals,statenames,'c50_per_state')
makepsychbarplot(rmax,animals,statenames,'rmax_per_state')
makepsychbarplot(baseline,animals,statenames,'baseline_per_state')
n=length(animals);
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

%% hyperbolic ratio fit
close all
[model_lowpupq,~,~]=HyFit([0 2 5 10 20 40 100].',M(:,1)); 
[model_highpupq,~,~]=HyFit([0 2 5 10 20 40 100].',M(:,2)); 
[model_highpupl,~,~]=HyFit([0 2 5 10 20 40 100].',M(:,3)); 

predicted_lowpupq=nanmean(predint(model_lowpupq,0:0.01:100),2); %predict values from function fit 
predicted_highpupq=nanmean(predint(model_highpupq,0:0.01:100),2);
predicted_highpupl=nanmean(predint(model_highpupl,0:0.01:100),2);
figure;
semilogx(0:0.01:100,predicted_lowpupq,'color','k'); %overlay early and late days psyc curve on semi log plot
hold on;
semilogx(0:0.01:100,predicted_highpupq,'color',[0.9290 0.6940 0.1250]);
hold on
semilogx(0:0.01:100,predicted_highpupl,'color','r');
hold on
errorbar([0 2 5 10 20 40 100].',M(:,1), S(:,1),'color','k','LineStyle','none');
hold on
errorbar([0 2 5 10 20 40 100].',M(:,2), S(:,2),'color',[0.9290 0.6940 0.1250],'LineStyle','none');
hold on
errorbar([0 2 5 10 20 40 100].',M(:,3), S(:,3),'color','r','LineStyle','none');
axis([0,100,0,100]);
set(gca,'XTickLabel', [0 0.1 1 10 100])
xlabel('Contrast');
ylabel('% Correct');
title(strcat('Psyc Curve Across States'));
mysave(gcf, fullfile(outputfiggolder,['psych_curve_hyperbolicfit_states' loosestr]));


n = length(animals); 
M = nanmean(suc_rate, 2)*100;
S = nanstd(suc_rate, [], 2)/sqrt(n-1)*100;
figure;barwitherr(S, M)
set(gca,'XTickLabel',statenames);
ylabel('% Success Rate');
mysave(gcf, fullfile(outputfiggolder,['success_rate_by_state' loosestr]));

mean_mean_across_groups1=M;
std_mean_across_groups1=S;
figure;
CondColors=[0,0,0;0.9290 0.6940 0.1250;1,0,0];
subplot(1,1,1)
set(gcf,'renderer','Painters')
hold on
for b = 1:3
    bg=bar(b, mean_mean_across_groups1(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
set(gca,'xtick',1:3)
set(gca,'xticklabel',statenames)
h = errorbar(1:3,mean_mean_across_groups1, std_mean_across_groups1,'LineStyle','none','LineWidth',0.5);title('State Success Rate');
h.Color='k';
set(h, 'marker', 'none'); 
mysave(gcf, fullfile(outputfiggolder,['colored_success_rate_by_state' loosestr]), 'all');
end


function plot_prediction(isloose, animals, statenames, outputfiggolder)

if isloose
    loosestr = 'loose';
else
    loosestr = '';
end
for ai = 1:length(animals)
    Lab = load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animals{ai},'\',animals{ai},'trials_3states' loosestr '.mat'],...
    'low_pup_q','high_pup_q','high_pup_l','days_to_process', 't'); %#ok<NASGU>
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
            
        
    y = Lab.(statenames{state_i}).trialslabels.blinksummary;
    inds = y<3;
    labels = y(inds);
    
    
    
            slope_correct.(statenames{state_i})(:,ai) = mean(slopeData(:, labels==1),2);
            slope_incorrect.(statenames{state_i})(:,ai) = mean(slopeData(:, labels==2),2);
        end
    end
end

parcels_names = get_allen_meta_parcels;
n = length(animals); 
CondColors=[0,0,0;0.9290 0.6940 0.1250;1,0,0];
for state_i = 1:length(statenames)
    subplot(3,1,state_i);set(gcf,'renderer','Painters');
    M(:,1, state_i) = mean(slope_correct.(statenames{state_i}),2);
    M(:,2, state_i) = mean(slope_incorrect.(statenames{state_i}),2);
    S(:,1, state_i) = std(slope_correct.(statenames{state_i}),[],2)/sqrt(n-1);
    S(:,2, state_i) = std(slope_incorrect.(statenames{state_i}),[],2)/sqrt(n-1);
    h=barwitherr(S(:,:, state_i),M(:,:, state_i));title(statenames{state_i});
    h(1).EdgeColor = 'none';
    h(2).EdgeColor = 'none';
    set(h(1),'FaceColor','g');
    set(h(2),'FaceColor','r'); 
    set(h(1),'FaceAlpha',0.4);
    set(h(2),'FaceAlpha',0.4); 
    set(gca,'xtick',1:23)
    set(gcf, 'Position',  [1,1, 700,1000]);    
    set(gca,'xticklabel',parcels_names)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);
end
mysave(gcf,fullfile(outputfiggolder, ['slope_by_state_by_behavior' loosestr]));

plot_bars_3colors(squeeze(M(:,1,:)), squeeze(S(:,1,:)), statenames, parcels_names)
mysave(gcf,fullfile(outputfiggolder, ['slope_by_state_correct' loosestr]));
plot_bars_3colors(squeeze(M(:,2,:)), squeeze(S(:,2,:)), statenames, parcels_names)
mysave(gcf,fullfile(outputfiggolder, ['slope_by_state_incorrect' loosestr]));
clear M;
clear S;

figure;
mean_mean_across_groups1=nanmean(acc_all_parcels,2);
std_mean_across_groups1=nanstd(acc_all_parcels,[],2)/sqrt(n-1);
figure;
CondColors=[0,0,0;0.9290 0.6940 0.1250;1,0,0];
subplot(1,1,1)
set(gcf,'renderer','Painters')
hold on
for b = 1:3
    bg=bar(b, mean_mean_across_groups1(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
ylim([0.5 0.9]);
set(gca,'xtick',1:3)
set(gca,'xticklabel',statenames)
h = errorbar(1:3,mean_mean_across_groups1, std_mean_across_groups1,'LineStyle','none','LineWidth',0.5);title('Acc All Parcels Per State');
h.Color='k';
set(h, 'marker', 'none'); 
mysave(gcf,fullfile(outputfiggolder, ['behavior_prediction_by_state_all_parcels' loosestr]));

figure;
set(gcf,'renderer','Painters')
M = nanmean(acc_per_parcel,3);
S = nanstd(acc_per_parcel,[],3)/sqrt(n-1);
h=barwitherr(S,M);
CondColors=[0,0,0;0.9290 0.6940 0.1250;1,0,0];
h(1).EdgeColor = 'none';
h(2).EdgeColor = 'none';
h(3).EdgeColor = 'none';
set(h(1),'FaceColor',CondColors(1,:));
set(h(2),'FaceColor',CondColors(2,:));
set(h(3),'FaceColor',CondColors(3,:));
set(h(1),'FaceAlpha',0.8);
set(h(2),'FaceAlpha',0.8);
set(h(3),'FaceAlpha',0.8);
set(gcf, 'Position',  [150,150, 1500,700]);
legend(statenames)
ylim([0.5 0.9]);
set(gca,'XTick', 1:length(parcels_names));
set(gca,'XTickLabel', parcels_names);

set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
set(gca,'XTickLabelRotation',45);
set(gcf,'Position',[1          41        1500         700])
mysave(gcf,fullfile(outputfiggolder, ['behavior_prediction_by_state_per_parcel' loosestr]));

%%difference plots
%pupil high (2)-low(1)
%run(3)-low(1)
%run(3)-pupilhigh(2)
braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');

graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'),M(:,2),M(:,1),'trials/spon_pupilhigh_pupillow/','accuracy_SVM_diff','pupil high - low (svm acc)',parcels_names,n)
graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
    M(:,2),M(:,1),'trials/spon_pupilhigh_pupillow/','accuracy_SVM_heatmap','pupil high - low (svm acc)');

graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'),M(:,3),M(:,1),'trials/spon_run_pupillow/','accuracy_SVM_diff','run - pupil low (svm acc)',parcels_names,n)
graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
    M(:,3),M(:,1),'trials/spon_run_pupillow/','accuracy_SVM_heatmap','run - pupil low (svm acc)');

graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'),M(:,3),M(:,2),'trials/spon_run_pupilhigh/','accuracy_SVM_diff','run - pupil high (svm acc)',parcels_names,n)
graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
    M(:,3),M(:,2),'trials/spon_run_pupilhigh/','accuracy_SVM_heatmap','run - pupil high (svm acc)');
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

function makepsychbarplot(c50s,animals,statenames,name)
n = length(animals); 
mean_mean_across_groups1=nanmean(c50s,1);
std_mean_across_groups1=nanstd(c50s,[],1)./sqrt(n-1);
figure;
CondColors=[0,0,0;0.9290 0.6940 0.1250;1,0,0];
subplot(1,1,1)
set(gcf,'renderer','Painters')
hold on
for b = 1:3
    bg=bar(b, mean_mean_across_groups1(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
set(gca,'xtick',1:3)
set(gca,'xticklabel',statenames)
h = errorbar(1:3,mean_mean_across_groups1, std_mean_across_groups1,'LineStyle','none','LineWidth',0.5);title(name);
h.Color='k';
set(h, 'marker', 'none'); 
mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\','figure_1',name), 'all');
end







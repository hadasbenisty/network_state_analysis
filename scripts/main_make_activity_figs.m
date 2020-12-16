function main_make_activity_figs
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));


daysmax=30;
animals={'xw','xx','xz','xt','xu' 'xs'};
stateslabels = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
outputfiggolder = 'X:\Hadas\Meso-imaging\lan\meso_results\figs\activity';
mkNewDir(outputfiggolder);
contrast_levels = [0 2 5 10 20 40 100];
spatialindex=[1 14 23];
%% Fig 2 - activity
plot_slopes(spatialindex, animals, stateslabels, outputfiggolder,daysmax )
makeslopeamplitudeplots(spatialindex, animals, stateslabels, outputfiggolder, contrast_levels,daysmax);
plotCRF(outputfiggolder, stateslabels, contrast_levels, animals,daysmax);
close all;
end


function makeslopeamplitudeplots(spatialindex, animals, statenames, outputfiggolder, contrast_levels,daysmax)
Ntime = 100;
tN = linspace(-4,4,Ntime);
x.corr = nan(length(spatialindex), Ntime,length(animals), length(contrast_levels));
x.incorr = nan(length(spatialindex),Ntime,length(animals), length(contrast_levels));
for state_i=1:length(statenames)
    
    for animal_i=1:length(animals)
        animal=char(animals(animal_i));
        load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_meta_data']),...
            'behavior_labels', 'arousal_labels', 'contrast_labels', 'days_labels');
        load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_imaging_time_traces_global_Allen_dfff.mat']),'t','imaging_time_traces');
        for ci=1:length(contrast_levels)
            indxcorr=behavior_labels==1& arousal_labels==state_i & ...
                contrast_labels== contrast_levels(ci)&days_labels<=daysmax ;
            if sum(indxcorr) > 0
            curr_corr=imaging_time_traces(spatialindex,:,indxcorr);
            vectorcorr=interp1(t, mean(curr_corr,3)', tN)';
            x.corr(:, :,animal_i, ci)=vectorcorr;
            end
            indxincorr=behavior_labels==2& arousal_labels==state_i & ...
                contrast_labels== contrast_levels(ci)&days_labels<=daysmax ;
            if sum(indxincorr) > 0
                curr_incorr=imaging_time_traces(spatialindex,:,indxincorr);
                vectorincorr=interp1(t, mean(curr_incorr,3)', tN)';
                x.incorr(:, :,animal_i, ci)=vectorincorr;
            end
        end
    end
    parcelnames = get_allen_meta_parcels;
    stind=findClosestDouble(tN,-0.5);enind=findClosestDouble(tN,2);
    %% plot correct and incorrect
    x1 = 0.0; x2 = 0.500;
    y1 = -2; y2 = 12;
    for i=1:length(spatialindex)
        for ci=1:length(contrast_levels)
            figure;
            set(gcf,'renderer','painters');
            fill([x1 x1 x2 x2],[y1 y2 y2 y1],[0.8 0.8 0.8],'LineStyle','none')
            hold on
            x3 = 0.450; x4 = 0.500;
            fill([x3 x3 x4 x4],[y1 y2 y2 y1],[0.3020 0.7490 0.9294],'LineStyle','none')
            xlim([-0.5 2]);ylim([-1 2])
            hold on
            shadedErrorBar(transpose(tN(stind:enind)),nanmean(x.corr((i), stind:enind,:,ci),3),(nanstd(x.corr((i),stind:enind,:,ci),[],3)./(sqrt(length(animals)-1))),'lineprops','g');
            hold on
            shadedErrorBar(transpose(tN(stind:enind)),nanmean(x.incorr((i),stind:enind,:,ci),3),(nanstd(x.incorr((i),stind:enind,:,ci),[],3)./(sqrt(length(animals)-1))),'lineprops','r');
            xlabel('Time [sec]');ylabel('Z-DF/F');title(strcat('Avg Correct vs Incorrect',statenames{state_i}));
            hold off
            mysave(gcf, fullfile(outputfiggolder,[parcelnames{spatialindex(i)},'_', num2str(contrast_levels(ci)),'contrast_corr_incorr_',statenames{state_i}]), 'all');
        end
    end
end
end

function plotCRF(outputfiggolder, statenames, contrasts, animals, daysmax)
acs=NaN(length(animals),length(statenames),length(contrasts));

for animal_i=1:length(animals)
    animal=animals{animal_i};
    load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_meta_data']),...
        'behavior_labels', 'contrast_labels', 'arousal_labels','days_labels');
    load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_imaging_time_traces_global_Allen_dfff.mat']),'t','imaging_time_traces');
    
    for state_i=1:length(statenames)
        %concatenate mean response per animal in 300 ms
        for contrast_i=1:length(contrasts)
            indx=behavior_labels<3& arousal_labels==state_i & ...
                contrast_labels== contrasts(contrast_i)&days_labels<=daysmax ;
            v1_act=imaging_time_traces(2,:,indx);
            %row wise zscoring
            st=findClosestDouble(t,0);
            ed=findClosestDouble(t,0.3);
            acs(animal_i,state_i,contrast_i)=nanmean(nanmean(v1_act(1,st:ed,:),2),3);
        end
    end
end
CondColors=get_3states_colors; 
figure;
xlim([0 100]);
hold on

for state_i=1:length(statenames)
    M=squeeze(acs(:,state_i,:));
    currmean=nanmean(M,1);%overlay early and late days psyc curve on semi log plo
    currsem=(nanstd(M,[],1)./(sqrt(length(animals)-1)));
    errorbar(contrasts,currmean,currsem,'color',CondColors(state_i,:),'LineWidth',1.5);
    CRFinfo.mean_per_state.(statenames{state_i})=currmean(2:7);
    CRFinfo.sem_per_state.(statenames{state_i})=currsem(2:7);
    hold on
    clearvars currmean currsem
end

CRFinfo.contrasts=contrasts(2:7);
New_XTickLabel = get(gca,'xtick');title('CRF per state');ylabel('Z-dF/F')
set(gca,'XTickLabel',New_XTickLabel);legend(statenames);xlabel('Contrast');
mysave(gcf, fullfile(outputfiggolder, 'CRF_Zscored_perstate'), 'all');
end





function plot_slopes(spatialindex, animals, statenames, outputfiggolder,daysmax)
mkNewDir(outputfiggolder);
slope_trial_time_start = 0.1;
slope_trial_time_end = 0.33;

if exist('X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\slopes.mat', 'file')
    
    load('X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\slopes.mat',...
        'slope_correct','slope_incorrect');
else
    for ai = 1:length(animals)
        load(['X:\Hadas\Meso-imaging\lan\' animals{ai} 'psych\spt\' animals{ai} '_trial_meta_data'],...
            'behavior_labels', 'arousal_labels','days_labels');
        
        
        load(['X:\Hadas\Meso-imaging\lan\' animals{ai} 'psych\spt\' animals{ai} '_trial_imaging_time_traces_global_Allen_dfff.mat'], ...
            'imaging_time_traces', 't' );
        
        for state_i = 1:length(statenames)
            
            
            
            inds = behavior_labels<3 &  arousal_labels==state_i & days_labels<=daysmax;
            labels = behavior_labels(inds);
            curr_imaging_time_traces=imaging_time_traces(:,:,inds);
            slopeData=nan(size(curr_imaging_time_traces,1), size(curr_imaging_time_traces,3));
            for T=1:size(curr_imaging_time_traces,3)
                slopeData(:,T) = getSlope(curr_imaging_time_traces(:, t>=slope_trial_time_start & t<slope_trial_time_end, T));
            end
            slope_correct.(statenames{state_i})(:,ai) = mean(slopeData(:, labels==1),2);
            slope_incorrect.(statenames{state_i})(:,ai) = mean(slopeData(:, labels==2),2);
        end
    end

    save('X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\slopes.mat',...
        'slope_correct','slope_incorrect');
end
parcels_names = get_allen_meta_parcels;
strstatenames=cell(3,1);
for k=1:length(statenames)
    strstatenames{k} =  statenames{k};
    strstatenames{k}(statenames{k}=='_') = ' ';
end



M=nan(length(parcels_names), 2, length(statenames));
S=nan(length(parcels_names), 2, length(statenames));
n = length(animals);
for k=1:length(statenames)
    M(:,1,k)=nanmean(slope_correct.(statenames{k}),2);
    M(:,2,k)=nanmean(slope_incorrect.(statenames{k}),2);
    S(:,1,k)=nanstd(slope_correct.(statenames{k}),[],2)/sqrt(n-1);
    S(:,2,k)=nanstd(slope_incorrect.(statenames{k}),[],2)/sqrt(n-1);
end
plot_correct_incorrect_per_3parcels(M, S, parcels_names, statenames,spatialindex);
legend('Correct','Incorrect')
mysave(gcf,fullfile(outputfiggolder, 'slope_per_state_correct_incorrect_3_parcels'));

end
function plot_correct_incorrect_per_3parcels(M, S, parcels_names, statenames,spatialindex)

CondColors = get_3states_colors;
for s=1:length(statenames)
    statenames{s}(statenames{s}=='_') = ' ';
end
figure;
set(gcf,'renderer','painters');
for parcel_i = 1:length(spatialindex)
    subplot(length(spatialindex),1,(parcel_i));
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');hold all;
    h1(1).YData(2:3)=0;
    h1(2).YData(2:3)=0;
    h1(1).FaceColor=CondColors(1,:);
    h1(2).FaceColor=CondColors(1,:);
    h1(2).FaceAlpha=0.4;
    h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');hold all;
    h1(1).YData([1 3])=0;
    h1(2).YData([1 3])=0;
    h1(1).FaceColor=CondColors(2,:);
    h1(2).FaceColor=CondColors(2,:);
    h1(2).FaceAlpha=0.4;h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');
    h1(1).YData([1 2])=0;
    h1(2).YData([1 2])=0;
    h1(1).FaceColor=CondColors(3,:);
    h1(2).FaceColor=CondColors(3,:);
    h1(2).FaceAlpha=0.4;h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    
    title(parcels_names{spatialindex(parcel_i)});
    
    set(gca,'xtick',1:23)
    set(gcf, 'Position',  [1,1, 700,1000]);
    set(gca,'xticklabel',statenames)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);%ylim([0 400]);
end
end




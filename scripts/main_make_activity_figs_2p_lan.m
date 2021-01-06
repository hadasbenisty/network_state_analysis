function main_make_activity_figs_2p_lan
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));


daysmax=30;
T = readtable('X:\Lan\FMB208 conditioning imaging data\animal_log_imaging_lt.xlsx');
animals = T.AnimalID;
poptype = T.TargetPop;
i = find(strcmp(animals, 'zb'));
i(2) = find(strcmp(animals, 'zy'));
i(3) = find(strcmp(animals, 'zn'));
i(4) = find(strcmp(animals, 'zi'));
i(5) = find(strcmp(animals, 'zc'));
i(6) = find(strcmp(animals, 'zm'));

animals=animals(setdiff(1:length(animals), i));
poptype=poptype(setdiff(1:length(poptype), i));
stateslabels = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
outputfiggolder = 'X:\Hadas\Meso-imaging\lan\meso_results\figs2p\activity';
mkNewDir(outputfiggolder);
contrast_levels = [0 2 5 10 20 40 100];
%% Fig 2 - activity
plot_slopes(animals, poptype, stateslabels, outputfiggolder,daysmax )
% makeslopeamplitudeplots(animals, poptype, stateslabels, outputfiggolder, contrast_levels,daysmax);
plotCRF(outputfiggolder, poptype, stateslabels, contrast_levels, animals,daysmax);
plotCRFbyslopes(outputfiggolder, poptype, stateslabels, contrast_levels, animals,daysmax)
close all;
end


function makeslopeamplitudeplots(animalsAll, poptype, statenames, outputfiggolder, contrast_levels,daysmax)
statenamesStr = {'Low Pupil Q', 'High Pupil Q','High Pupil L'};
Ntime = 100;
tN = linspace(-4,4,Ntime);
populations = unique(poptype);
for population_i = 1:length(populations)
    currpop = populations{population_i};
    animals = animalsAll(strcmp(poptype, currpop));
    
    for state_i=1:length(statenames)
        x.corr = nan(Ntime,length(animals), length(contrast_levels));
        x.incorr = nan(Ntime,length(animals), length(contrast_levels));
        x.corrAllc = nan(Ntime,length(animals));
        x.incorrAllc = nan(Ntime,length(animals));
        for animal_i=1:length(animals)
            animal=char(animals(animal_i));
            load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_meta_data']),...
                'behavior_labels', 'arousal_labels', 'contrast_labels', 'days_labels');
            load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_imaging_time_traces_global_dfff.mat']),'t','imaging_time_tracesN');
            
            
            indxcorr=behavior_labels==1& arousal_labels==state_i & ...
                days_labels<=daysmax ;
            if sum(indxcorr) > 0
                curr_corr=mean(imaging_time_tracesN(:,:,indxcorr));
                vectorcorr=interp1(t, mean(curr_corr,3)', tN)';
                x.corrAllc(:, animal_i)=vectorcorr;
            end
            indxincorr=behavior_labels==2& arousal_labels==state_i & ...
                days_labels<=daysmax ;
            if sum(indxincorr) > 0
                curr_incorr=mean(imaging_time_tracesN(:,:,indxincorr));
                vectorincorr=interp1(t, mean(curr_incorr,3)', tN)';
                x.incorrAllc(:, animal_i)=vectorincorr;
            end
            
            for ci=1:length(contrast_levels)
                indxcorr=behavior_labels==1& arousal_labels==state_i & ...
                    contrast_labels== contrast_levels(ci)&days_labels<=daysmax ;
                if sum(indxcorr) > 0
                    curr_corr=mean(imaging_time_tracesN(:,:,indxcorr));
                    vectorcorr=interp1(t, mean(curr_corr,3)', tN)';
                    x.corr(:, animal_i, ci)=vectorcorr;
                end
                indxincorr=behavior_labels==2& arousal_labels==state_i & ...
                    contrast_labels== contrast_levels(ci)&days_labels<=daysmax ;
                if sum(indxincorr) > 0
                    curr_incorr=mean(imaging_time_tracesN(:,:,indxincorr));
                    vectorincorr=interp1(t, mean(curr_incorr,3)', tN)';
                    x.incorr(:, animal_i, ci)=vectorincorr;
                end
            end
        end
        stind=findClosestDouble(tN,-3);enind=findClosestDouble(tN,2);
        %% plot correct and incorrect
        x1 = 0.0; x2 = 0.500;
        y1 = -1; y2 = 12;
        figure;
        for ci=1:length(contrast_levels)
            
            subplot(2,4,ci);
            set(gcf,'renderer','painters');
            fill([x1 x1 x2 x2],[y1 y2 y2 y1],[0.8 0.8 0.8],'LineStyle','none')
            hold on
            x3 = 0.450; x4 = 0.500;
            fill([x3 x3 x4 x4],[y1 y2 y2 y1],[0.3020 0.7490 0.9294],'LineStyle','none')
            xlim([-1 2]);ylim([-1 3])
            hold on
            shadedErrorBar(transpose(tN(stind:enind)),nanmean(x.corr(stind:enind,:,ci),2),(nanstd(x.corr(stind:enind,:,ci),[],2)./(sqrt(length(animals)-1))),'lineprops','g');
            hold on
            shadedErrorBar(transpose(tN(stind:enind)),nanmean(x.incorr(stind:enind,:,ci),2),(nanstd(x.incorr(stind:enind,:,ci),[],2)./(sqrt(length(animals)-1))),'lineprops','r');
            xlabel('Time [sec]');ylabel('Z-score');
            title([' C=', num2str(contrast_levels(ci))]);
            hold off
        end
        subplot(2,4,8);
        set(gcf,'renderer','painters');
        fill([x1 x1 x2 x2],[y1 y2 y2 y1],[0.8 0.8 0.8],'LineStyle','none')
        hold on
        x3 = 0.450; x4 = 0.500;
        fill([x3 x3 x4 x4],[y1 y2 y2 y1],[0.3020 0.7490 0.9294],'LineStyle','none')
        xlim([-1 2]);ylim([-1 3])
        hold on
        shadedErrorBar(transpose(tN(stind:enind)),nanmean(x.corrAllc(stind:enind,:),2),(nanstd(x.corr(stind:enind,:),[],2)./(sqrt(length(animals)-1))),'lineprops','g');
        hold on
        shadedErrorBar(transpose(tN(stind:enind)),nanmean(x.incorr(stind:enind,:),2),(nanstd(x.incorr(stind:enind,:),[],2)./(sqrt(length(animals)-1))),'lineprops','r');
        xlabel('Time [sec]');ylabel('Z-score');
        title(['All Contrasts']);
        hold off
        
        suptitle([currpop ' Correct vs Incorrect on ',statenamesStr{state_i}]);
        mysave(gcf, fullfile(outputfiggolder,[currpop,'_','contrast_corr_incorr_',statenames{state_i}]), 'all');
        
        
    end
end
end
function plotCRFbyslopes(outputfiggolder, poptype, statenames, contrasts, animalsAll,daysmax)
populations = unique(poptype);
for population_i = 1:length(populations)
    currpop = populations{population_i};
    animals = animalsAll(strcmp(poptype, currpop));
    acs=NaN(length(animals),length(statenames),length(contrasts));
    
    for animal_i=1:length(animals)
        animal=animals{animal_i};
        load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_meta_data']),...
            'behavior_labels', 'contrast_labels', 'arousal_labels','days_labels');
        load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_imaging_time_traces_global_dfff.mat']),'t','imaging_time_tracesN');
       
        for state_i=1:length(statenames)
            %concatenate mean response per animal in 300 ms
            for contrast_i=1:length(contrasts)
                indx=behavior_labels<3& arousal_labels==state_i & ...
                    contrast_labels== contrasts(contrast_i)&days_labels<=daysmax ;
                curr_imaging_time_traces=imaging_time_tracesN(:,:,indx);
                slopeData = nan(size(curr_imaging_time_traces,1), size(curr_imaging_time_traces,3));
                for T=1:size(curr_imaging_time_traces,3)
                    slopeData(:,T) = getSlope(curr_imaging_time_traces(:, t>=0 & t<0.33, T));
                end
                
                acs(animal_i,state_i,contrast_i)=nanmean(slopeData(:));
            end
        end
    end
    %     acs = acs - acsbase;
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
    New_XTickLabel = get(gca,'xtick');title('CRF per state');ylabel('Z-score')
    set(gca,'XTickLabel',New_XTickLabel);legend(statenames);xlabel('Contrast');
    mysave(gcf, fullfile(outputfiggolder, ['CRF_Zscored_perstate' currpop 'byslopes']), 'all');
end

end

function plotCRF(outputfiggolder, poptype, statenames, contrasts, animalsAll,daysmax)
populations = unique(poptype);
for population_i = 1:length(populations)
    currpop = populations{population_i};
    animals = animalsAll(strcmp(poptype, currpop));
    acs=NaN(length(animals),length(statenames),length(contrasts));
    
    for animal_i=1:length(animals)
        animal=animals{animal_i};
        load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_meta_data']),...
            'behavior_labels', 'contrast_labels', 'arousal_labels','days_labels');
        load(fullfile(['X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'],[animal '_trial_imaging_time_traces_global_dfff.mat']),'t','imaging_time_tracesN');
        
        for state_i=1:length(statenames)
            %concatenate mean response per animal in 300 ms
            for contrast_i=1:length(contrasts)
                indx=behavior_labels<3& arousal_labels==state_i & ...
                    contrast_labels== contrasts(contrast_i)&days_labels<=daysmax ;
                
                
                
                v1_act=mean(imaging_time_tracesN(:,:,indx),1);
                %row wise zscoring
                st=findClosestDouble(t,0);
                ed=findClosestDouble(t,0.3);
                acs(animal_i,state_i,contrast_i)=nanmean(nanmean(v1_act(1,st:ed,:),2),3);
                st=findClosestDouble(t,-.3);
                ed=findClosestDouble(t,0);
                acsbase(animal_i,state_i,contrast_i)=nanmean(nanmean(v1_act(1,st:ed,:),2),3);
            end
        end
    end
    %     acs = acs - acsbase;
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
    New_XTickLabel = get(gca,'xtick');title('CRF per state');ylabel('Z-score')
    set(gca,'XTickLabel',New_XTickLabel);legend(statenames);xlabel('Contrast');
    mysave(gcf, fullfile(outputfiggolder, ['CRF_Zscored_perstate' currpop]), 'all');
end

end



function plot_slopes(animalsAll, poptype, statenames, outputfiggolder,daysmax)
mkNewDir(outputfiggolder);
slope_trial_time_start = 0.1;
slope_trial_time_end = 0.33;
populations = unique(poptype);
for population_i = 1:length(populations)
    currpop = populations{population_i};
    animals = animalsAll(strcmp(poptype, currpop));
    if exist(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\slopes' currpop '.mat'], 'file')
        
        load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\slopes' currpop '.mat'],...
            'slope_correct','slope_incorrect');
    else
        
        for ai = 1:length(animals)
            load(['X:\Hadas\Meso-imaging\lan\' animals{ai} 'psych\spt\' animals{ai} '_trial_meta_data'],...
                'behavior_labels', 'arousal_labels','days_labels');
            
            
            load(['X:\Hadas\Meso-imaging\lan\' animals{ai} 'psych\spt\' animals{ai} '_trial_imaging_time_traces_global_dfff.mat'], ...
                'imaging_time_tracesN', 't' );
            
            
            for state_i = 1:length(statenames)
                
                
                
                inds = behavior_labels<3 &  arousal_labels==state_i & days_labels<=daysmax;
                labels = behavior_labels(inds);
                curr_imaging_time_traces=imaging_time_tracesN(:,:,inds);
                slopeData=nan(size(curr_imaging_time_traces,1), size(curr_imaging_time_traces,3));
                for T=1:size(curr_imaging_time_traces,3)
                    slopeData(:,T) = getSlope(curr_imaging_time_traces(:, t>=slope_trial_time_start & t<slope_trial_time_end, T));
                end
                
                slope_correct.(statenames{state_i})(1,ai) = mean(mean(slopeData(:, labels==1),2),1);
                slope_incorrect.(statenames{state_i})(1,ai) = mean(mean(slopeData(:, labels==2),2),1);
            end
        end
        
        save(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\slopes' currpop '.mat'],...
            'slope_correct','slope_incorrect');
    end
end
strstatenames=cell(3,1);
for k=1:length(statenames)
    strstatenames{k} =  statenames{k};
    strstatenames{k}(statenames{k}=='_') = ' ';
end
M=nan(length(populations), 2, length(statenames));
S=nan(length(populations), 2, length(statenames));
for population_i = 1:length(populations)
    currpop = populations{population_i};
    load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\slopes' currpop '.mat'],...
        'slope_correct','slope_incorrect');
    
    
    
    
    
    n = length(animals);
    for k=1:length(statenames)
        M(population_i,1,k)=nanmean(slope_correct.(statenames{k}),2);
        M(population_i, 2,k)=nanmean(slope_incorrect.(statenames{k}),2);
        S(population_i, 1,k)=nanstd(slope_correct.(statenames{k}),[],2)/sqrt(n-1);
        S(population_i, 2,k)=nanstd(slope_incorrect.(statenames{k}),[],2)/sqrt(n-1);
    end
end

h=plot_correct_incorrect_per_state_per_parcels(M, S, populations, strstatenames);
for k=1:length(h)
    ylabel(h(k), 'Slope');xlabel(h(k), 'Cell Population'); legend(h(k), 'Correct','Incorrect')
end

mysave(gcf,fullfile(outputfiggolder, ['slope_per_state_correct_incorrect' ]));



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




% Fig 1 - Behavior
function main_make_behavior_figs_2p_lan

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
outputfiggolder = 'X:\Hadas\Meso-imaging\lan\meso_results\figs2p\behavior';
contrast_levels = [0 2 5 10 20 40 100];
daysmax=30;
% plot_time_spent(animals, poptype, stateslabels, outputfiggolder,daysmax)
plot_psych_curve_per_state(animals, poptype, stateslabels, contrast_levels, outputfiggolder,daysmax)

close all

end






function plot_psych_curve_per_state(animalsAll, poptype, statenames, contrast_levels, outputfiggolder,daysmax)
populations = unique(poptype);
for population_i = 1:length(populations)
    currpop = populations{population_i};
    animals = animalsAll(strcmp(poptype, currpop));
suc_rate=nan(length(statenames), length(animals));
if exist(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\psych_curv_data' currpop '.mat'], 'file')
    load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\psych_curv_data' currpop '.mat'], 'psych_curv', 'suc_rate');
else
    psych_curv=nan(length(contrast_levels), length(statenames), length(animals));
    for ai=1:length(animals)
        if ~exist(fullfile('X:\Hadas\Meso-imaging\lan\', [animals{ai} 'psych'], 'spt',[animals{ai},'_trial_meta_data.mat']), 'file')
            continue;
        end
        load(fullfile('X:\Hadas\Meso-imaging\lan\', [animals{ai} 'psych'], 'spt',[animals{ai},'_trial_meta_data']),...
            'behavior_labels','contrast_labels','arousal_labels','days_labels');
        
        
        for state_i = 1:length(statenames)
            disp(statenames{state_i})
            
            
            inds = behavior_labels<3&arousal_labels==state_i&days_labels<=daysmax;
            labels = behavior_labels(inds);
            contrast = contrast_labels(inds);
            suc_rate(state_i, ai) = sum(labels==1)/length(labels);
            for ci=1:length(contrast_levels)
                psych_curv(ci, state_i, ai) =  sum(labels==1 & contrast == contrast_levels(ci))/sum(contrast == contrast_levels(ci));
            end
        end
    end
    save(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\psych_curv_data' currpop '.mat'], 'psych_curv', 'suc_rate');
    
end
%% c50 and rmax per animal for barplots

c50s=NaN(length(animals),length(statenames));
rmax=NaN(length(animals),length(statenames));
baseline=NaN(length(animals),length(statenames));
discardnans = isnan(squeeze(sum(sum(psych_curv))));
for animal_i=1:length(animals)
    for state_i=1:length(statenames)
        [~,coefs]=HyFit([0 2 5 10 20 40 100].',psych_curv(:,state_i,animal_i));
        c50s(animal_i,state_i)=coefs(3);
        baseline(animal_i,state_i)=coefs(4);
        rmax(animal_i,state_i)=coefs(1);
        clearvars coefs
    end
end
statenamesstr = statenames;
for si=1:length(statenames)
    statenamesstr{si}(statenames{si}=='_') = ' ';
end
mkNewDir(outputfiggolder);
makepsychbarplot(c50s,animals,statenamesstr,'c50')
mysave(gcf, fullfile(outputfiggolder,['c50_per_state' currpop]));
makepsychbarplot(rmax,animals,statenamesstr,'R_{max}')
mysave(gcf, fullfile(outputfiggolder,['rmax_per_state' currpop]));
makepsychbarplot(baseline,animals,statenamesstr,'baseline')
mysave(gcf, fullfile(outputfiggolder,['baseline_per_state' currpop]));

n=length(animals);
M = nanmean(psych_curv, 3)*100;
S = nanstd(psych_curv, [], 3)/sqrt(n-1)*100;
CondColors = get_3states_colors;

figure;
for si=1:length(statenames)
    errorbar(contrast_levels, M(:,si), S(:,si), 'Color',CondColors(si,:));
    hold all;
end
xlabel('% Contrast');
ylabel('% Success Rate');
legend(statenamesstr);
mysave(gcf, fullfile(outputfiggolder,['psych_curves_by_state' currpop]));

%% hyperbolic ratio fit
close all
figure;
for si=1:length(statenames)
    model_perstate=HyFit(contrast_levels.',M(:,si));
    predicted=nanmean(predint(model_perstate,0:0.01:100),2); %predict values from function fit
    semilogx(0:0.01:100,predicted,'color',CondColors(si,:)); %overlay early and late days psyc curve on semi log plot
    hold on;
end
for si=1:length(statenames)
    errorbar(contrast_levels.',M(:,si), S(:,si),'color',CondColors(si,:),'LineStyle','none');
    
end
legend(statenamesstr)
axis([0,100,0,100]);
set(gca,'XTickLabel', [0 0.1 1 10 100])
xlabel('Contrast');
ylabel('% Correct');
title(strcat('Psyc Curve Across States'));
mysave(gcf, fullfile(outputfiggolder,['psych_curve_hyperbolicfit_states' currpop]));


n = length(animals);
M = nanmean(suc_rate, 2)*100;
S = nanstd(suc_rate, [], 2)/sqrt(n-1)*100;
figure;
plot_3_bars(M, S, statenamesstr);
ylim([0 70]);
set(gca,'XTickLabel',statenamesstr);
ylabel('% Success Rate');
mysave(gcf, fullfile(outputfiggolder, ['success_rate_by_state' currpop] ));
end
end




function makepsychbarplot(c50s,animals,statenames,name)
n = length(animals);
mean_mean_across_groups1=nanmean(c50s,1);
std_mean_across_groups1=nanstd(c50s,[],1)./sqrt(n-1);
figure;
CondColors = get_3states_colors;
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
% mysave(gcf, fullfile('X:\Lav\ProcessingDirectory2p\','figure_1',name), 'all');
end
function plot_time_spent(animalsAll, poptype, stateslabels, outputfiggolder,daysmax)
mkNewDir(outputfiggolder);
behavelabel = {'correct','incorrect'};
populations = unique(poptype);
for population_i = 1:length(populations)
    currpop = populations{population_i};
    animals = animalsAll(strcmp(poptype, currpop));
if exist(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\time_spent_data' currpop '.mat'], 'file')
    load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\time_spent_data' currpop '.mat'],'Nall',...
        'Ncor','Ninc','pupvals','wheelvals');
else
    for state_i = 1:length(stateslabels)
        for bi = 1:length(behavelabel)
            wheelvals.(stateslabels{state_i}).(behavelabel{bi})=[];
            pupvals.(stateslabels{state_i}).(behavelabel{bi})=[];
        end
    end
    
    Nall = nan(length(animals), length(stateslabels));
    Ncor = nan(length(animals), length(stateslabels));
    Ninc = nan(length(animals), length(stateslabels));
    
    for ai = 1:length(animals)
        if ~exist(fullfile('X:\Hadas\Meso-imaging\lan\', [animals{ai} 'psych'], 'spt',[animals{ai},'_trial_meta_data.mat']), 'file')
            continue;
        end
        data = load(fullfile('X:\Hadas\Meso-imaging\lan\', [animals{ai} 'psych'], 'spt',[animals{ai},'_trial_meta_data']));
        
        for state_i = 1:length(stateslabels)
            Nall(ai, state_i) = sum(data.behavior_labels<3&data.arousal_labels==state_i&data.days_labels<=daysmax);
            Ncor(ai, state_i) = sum(data.behavior_labels==1&data.arousal_labels==state_i&data.days_labels<=daysmax);
            Ninc(ai, state_i) = sum(data.behavior_labels==2&data.arousal_labels==state_i&data.days_labels<=daysmax);
            for bi = 1:length(behavelabel)
                x = data.wheel_trace(:,data.behavior_labels==bi&data.arousal_labels==state_i);
                wheelvals.(stateslabels{state_i}).(behavelabel{bi}) = cat(1, wheelvals.(stateslabels{state_i}).(behavelabel{bi}), x(:));
                x = data.pupil_trace(:,data.behavior_labels==bi&data.arousal_labels==state_i);
                
                pupvals.(stateslabels{state_i}).(behavelabel{bi}) = cat(1, pupvals.(stateslabels{state_i}).(behavelabel{bi}), x(:));
            end
        end
    end
    save(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\time_spent_data' currpop '.mat'],'Nall',...
        'Ncor','Ninc','pupvals','wheelvals');
end
n = sum(~isnan(sum(Nall,2)));
n=length(animals);
CondColors = get_3states_colors;
stateslabels_ttls = {'low pup q' 'high pup q' 'high pup l'};

%% wheel histograms per state and behavior
binsN = linspace(0,60,32);
figure;
l=1;
for state_i = 1:length(stateslabels)
    for bi = 1:length(behavelabel)
        subplot(length(stateslabels),length(behavelabel),l);
        set(gcf,'renderer','Painters');
        histogram(wheelvals.(stateslabels{state_i}).(behavelabel{bi}), binsN,'facecolor',CondColors(state_i,:),'facealpha',.8,'edgecolor','none');
        title([stateslabels_ttls{state_i} ' ' behavelabel{bi}]);xlim([0 60]); xlabel('wheel');
        l=l+1;
    end
end
mysave(gcf, fullfile(outputfiggolder, ['trial_wheel_hist_3states' currpop]));

%% pupil histograms per state and behavior
binsN = linspace(0,6e3,32);
figure;
l=1;
for state_i = 1:length(stateslabels)
    for bi = 1:length(behavelabel)
        subplot(length(stateslabels),length(behavelabel),l);
        set(gcf,'renderer','Painters');
        histogram(pupvals.(stateslabels{state_i}).(behavelabel{bi}), binsN,'facecolor',CondColors(state_i,:),'facealpha',.8,'edgecolor','none');
        title([stateslabels_ttls{state_i} ' ' behavelabel{bi}]);xlim([0 6e3]); xlabel('pupil');
        l=l+1;
    end
end
mysave(gcf, fullfile(outputfiggolder, ['trial_pup_histogram_3states' currpop]));

%% time spent
Mall = nanmean(Nall);
Sall = nanstd(Nall)/sqrt(n-1);
Mcor = nanmean(Ncor);
Scor = nanstd(Ncor)/sqrt(n-1);
Minc = nanmean(Ninc);
Sinc = nanstd(Ninc)/sqrt(n-1);
figure;subplot(3,1,1);
plot_3_bars(Mall, Sall, stateslabels_ttls)
ylim([0 500]);title('All Trials');
subplot(3,1,2);plot_3_bars(Mcor, Scor, stateslabels_ttls)
ylim([0 200]);title('Correct Trials');
subplot(3,1,3);
plot_3_bars(Minc, Sinc, stateslabels_ttls)
ylim([0 300]);
title('Incorrect Trials');
mysave(gcf, fullfile(outputfiggolder, ['time_spent_bytrials_absolute_numbers' currpop]));



n=length(animals);
Mall = nanmean(Nall./sum(Nall,2));
Sall = nanstd(Nall./sum(Nall,2))/sqrt(n-1);
Mcor = nanmean(Ncor./sum(Ncor,2));
Scor = nanstd(Ncor./sum(Ncor,2))/sqrt(n-1);
Minc = nanmean(Ninc./sum(Ninc,2));
Sinc = nanstd(Ninc./sum(Ninc,2))/sqrt(n-1);

figure;subplot(4,1,1);
plot_3_bars(Mall, Sall, stateslabels_ttls)
ylim([0 .7]);title('All Trials');


subplot(4,1,2);
plot_3_bars(Mcor, Scor, stateslabels_ttls);ylim([0 .7]);
title('Correct Trials');
subplot(4,1,3);plot_3_bars(Minc, Sinc, stateslabels)
ylim([0 .7]);
title('Incorrect Trials');
%%
n=length(animals);
Nall = Ncor + Ninc;
Mcor = nanmean(Ncor./Nall);
Scor = nanstd(Ncor./Nall)/sqrt(n-1);
Minc = nanmean(Ninc./Nall);
Sinc = nanstd(Ninc./Nall)/sqrt(n-1);

subplot(4,1,4);
barwitherr([Scor; Sinc]' , [Mcor;Minc]');ylim([0 .7]);
set(gca,'XTickLabel', stateslabels_ttls);title('Time Spent By State');legend('Correct', 'Incorrect');


mysave(gcf, fullfile(outputfiggolder, ['time_spent_bytrials' currpop]));

end

end


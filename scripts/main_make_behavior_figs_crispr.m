% Fig 1 - Behavior
function main_make_behavior_figs_crispr
addpath(genpath('..\meta_data_processing'));
animals_db = get_animals_meta_data_by_csv;
addpath(genpath('../functions'))
addpath(genpath('../utils'))
stateslabels3 = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
stateslabels2 = { 'qui', 'loc'};
stateslabels5 = { 'sit','loc','low_pupil','high_pupil','low_face','high_face'};
procfolder = 'X:\Hadas\Meso-imaging\CRISPR\analysis_results\';
outputfiggolder = 'X:\Hadas\Meso-imaging\CRISPR\Figures\arousal';
mkNewDir(outputfiggolder);
mkNewDir(procfolder);
% plot_arousal_by_onset_noppl(procfolder, animals_db,spike2path, outputfiggolder)
% plot_arousal_by_onset_withppl(procfolder, animals_db, spike2path, outputfiggolder);
plot_time_spent(animals_db, stateslabels5, procfolder, outputfiggolder)

plot_time_spent(animals_db, stateslabels3, procfolder, outputfiggolder)
plot_time_spent(animals_db, stateslabels2, procfolder, outputfiggolder)


end

function plot_arousal_by_onset_noppl(procfolder, animals_db,spike2path, outputfiggolder)
n = length(animals_db.animal_list);
for ai = 1:n
  
    if  animals_db.isgoodpupil_list(ai)~=find(strcmp(animals_db.isgoodpupil_lut,'GOOD'))||~isfile(fullfile(procfolder, animals_db.folder_list{ai}, ['arousal_2state_ITI_segemts.mat']))
        continue;
    end
disp(ai);
    load(fullfile(spike2path, animals_db.folder_list{ai},'smrx_signals_v4.mat'),'timing','channels_data');
    wheelN = interp1([1:length(channels_data.wheelspeed)]/5e3,channels_data.wheelspeed, timing.mesostart);
    wheelN=zscore(wheelN);
    wheelonsets = timing.wheelOn;wheeltime = timing.mesostart;
    bef=10;af = 20;fsample = round(1/median(diff(timing.mesostart)));
   Wx=[];
    for k=1:length(wheelonsets)
        if wheelonsets(k)-bef<wheeltime(1)
            continue;
        end
        if wheelonsets(k)+af>wheeltime(end)
            break;
        end
        i = findClosestDouble(wheeltime, wheelonsets(k));
        Wx(:,k) = wheelN(i-bef*fsample:i+af*fsample);
    end
    t = linspace(-bef,af,size(Wx,1));
    figure;
   subplot(2,2,1);imagesc(t,1:size(Wx,2),(Wx'));title('Wheel');xlabel('Time [sec]');colorbar;
   subplot(2,2,2); 
    plot(t,mean(Wx,2));axis tight;xlabel('Time [sec]');legend('W');

    str=animals_db.folder_list{ai};
    str(str == '/') = '_';str(str == '\') = '_';
    ttl = animals_db.animal_lut{animals_db.animal_list(ai)};
    ttl(ttl == '_') = ' ';
     subplot(2,1,2);
plot(wheeltime,wheelN);hold all;xlabel('Time [sec]');legend('W');axis tight;
for k=1:length(wheelonsets)
    line(wheelonsets(k)*[1 1],get(gca,'YLim'),'Color','k');
end
set_large_fig;legend('W','Wheel Onset');axis tight;colormap('jet');
    suptitle(ttl)
   mysave(gcf,fullfile(outputfiggolder,['arousal_traces_' str]));
end
end

function plot_arousal_by_onset_withppl(procfolder, animals_db,spike2path, outputfiggolder)
n = length(animals_db.animal_list);
for ai = 1:n
    if animals_db.isgoodpupil_list(ai)~=find(strcmp(animals_db.isgoodpupil_lut,'GOOD'))||~isfile(fullfile(procfolder, animals_db.folder_list{ai}, ['arousal_3state_traces.mat']))
        continue;
    end

    disp(ai);
    load(fullfile(procfolder, animals_db.folder_list{ai}, ['arousal_3state_traces.mat']),...
        'pupil_Norm', 'pupil_time', 'wheelN');%
    load(fullfile(spike2path, animals_db.folder_list{ai},'smrx_signals_v4.mat'),'timing');
    wheelonsets = timing.wheelOn;
    bef=10;af = 20;fsample = round(1/median(diff(pupil_time)));
    Px = [];Wx=[];
    for k=1:length(wheelonsets)
        if wheelonsets(k)-bef<pupil_time(1)
            continue;
        end
        if wheelonsets(k)+af>pupil_time(end)
            break;
        end
        i = findClosestDouble(pupil_time, wheelonsets(k));
        if i-bef*fsample < 1
            continue;
        end
        Px(:,k) = pupil_Norm(i-bef*fsample:i+af*fsample);
        Wx(:,k) = wheelN(i-bef*fsample:i+af*fsample);
    end
    t = linspace(-bef,af,size(Px,1));
    figure;
    subplot(3,2,1);imagesc(t,1:size(Px,2),(Px'));title('Pupil');xlabel('Time [sec]');colorbar;
   subplot(3,2,2);imagesc(t,1:size(Wx,2),(Wx'));title('Wheel');xlabel('Time [sec]');colorbar;
   subplot(3,2,3); 
   plot(t,mean(Px,2));hold all;
    plot(t,mean(Wx,2));axis tight;xlabel('Time [sec]');legend('P','W');
    subplot(3,2,4);
    plotEmbeddingWithColors([mean(Px,2) mean(Wx,2)],t);
    str=animals_db.folder_list{ai};
    str(str == '/') = '_';str(str == '\') = '_';
    ttl = animals_db.animal_lut{animals_db.animal_list(ai)};
    ttl(ttl == '_') = ' ';
    subplot(3,1,3);plot(pupil_time,pupil_Norm)
hold all
plot(pupil_time,wheelN);hold all;xlabel('Time [sec]');legend('P','W');axis tight;
for k=1:length(wheelonsets)
    line(wheelonsets(k)*[1 1],get(gca,'YLim'),'Color','k');
end
set_large_fig;legend('P','W','Wheel Onset');axis tight;
    suptitle(ttl)
   mysave(gcf,fullfile(outputfiggolder,['arousal_traces_' str]));
end
end

function plot_time_spent(animals_db, stateslabels, procfolder, outputfiggolder)


if exist(fullfile(procfolder, ['time_spent_data' num2str(length(stateslabels)) 'states.mat']), 'file')
    load(fullfile(procfolder, ['time_spent_data' num2str(length(stateslabels)) 'states.mat']), 'pupvalsall','wheelvalsall', 'Nstates','Tstates');
else
    n = length(animals_db.animal_list);
    for state_i = 1:length(stateslabels)
        wheelvalsall.(stateslabels{state_i})=nan(n,1);
        pupvalsall.(stateslabels{state_i})=nan(n,1);
        
    end
    Nstates = nan(n, length(stateslabels));
    Tstates = nan(n, length(stateslabels));
    for ai = 1:n
        disp(ai/n);
        if ~isfile(fullfile(procfolder, animals_db.folder_list{ai}, ['arousal_traces_states.mat']))
            continue;
        end
        pupvals=[];
         
        load(fullfile(procfolder, animals_db.folder_list{ai}, ['arousal_traces_states.mat']));%
        for state_i = 1:length(stateslabels)
            if isfield(segments_arousals, stateslabels{state_i})
                if isempty(segments_arousals.(stateslabels{state_i}))
                    Tstates(ai, state_i) = 0;
                    Nstates(ai, state_i)=0;
                else
                    
                    Tstates(ai, state_i) = sum(diff(segments_arousals.(stateslabels{state_i}),1,2));
                    Nstates(ai, state_i) = size(segments_arousals.(stateslabels{state_i}),1);
                end
%                 wheelvalsall.(stateslabels{state_i})(ai) = mean(wheelvals.(stateslabels{state_i}));
%                 if ~isempty(pupvals)
%                     pupvalsall.(stateslabels{state_i})(ai) = mean(pupvals.(stateslabels{state_i}));
%                 end
            end
        end
    end
    save(fullfile(procfolder, ['time_spent_data' num2str(length(stateslabels)) 'states.mat']),'Tstates','Nstates', 'wheelvalsall', 'pupvalsall');
end
stateslabels_ttls = stateslabels;
for si=1:length(stateslabels_ttls)
    stateslabels_ttls{si}(strfind(stateslabels_ttls{si},'_')) = ' ';
end

CondColors = get_3states_colors;
if length(stateslabels_ttls)==2
    CondColors=CondColors([1 3],:);
end
figure;
for ci = 1:length(animals_db.type_lut)
b(:,ci) = hist(animals_db.arousal_cluster_list(animals_db.type_list==ci&animals_db.arousal_cluster_list>0),1:3);
% b(:,ci)=b(:,ci)/sum(b(:,ci));
end
bar(b);set(gca,'XTickLabel', animals_db.arousal_cluster_lut(3:end));
legend(animals_db.type_lut);

%% time spent
Mbytype=[];Sbytype=[];nbytype=[];
for ci = 2:length(animals_db.arousal_cluster_lut)
    X=Tstates(animals_db.arousal_cluster_list==ci-1 ,:);
    X=bsxfun(@rdivide, X, nansum(X,2));
    Mbytype(ci,:) = nanmean(X);
    Sbytype(ci,:) = nanstd(X);
    nbytype(ci,:) = sum(~isnan(X));
    
end
figure;b = barwitherr(Sbytype./sqrt(nbytype-1), Mbytype);
xlim([1.5 4.5])
set(gca,'XTickLabel', animals_db.arousal_cluster_lut(2:end));
legend(stateslabels_ttls);
title('Fraction of time on arousal state');

Mbytype=[];Sbytype=[];nbytype=[];

for ci = 1:length(animals_db.type_lut)
    X=Tstates(animals_db.type_list==ci &animals_db.arousal_cluster_list>0, :);
    X=bsxfun(@rdivide, X, nansum(X,2));
    Mbytype(ci,:) = nanmean(X);
    Sbytype(ci,:) = nanstd(X);
    nbytype(ci,:) = sum(~isnan(X));
    
end
figure;b = barwitherr(Sbytype./sqrt(nbytype-1), Mbytype);
for i=1:length(stateslabels_ttls)
    b(i).FaceColor = CondColors(i,:);
end
set(gca,'XTickLabel', animals_db.type_lut);
legend(stateslabels_ttls);
title('Fraction of time on arousal state');

mysave(gcf, fullfile(outputfiggolder, ['time_spent_fraction', num2str(length(stateslabels_ttls)) ]));

figure;b = barwitherr(Sbytype'./sqrt(nbytype'-1), Mbytype');

legend( animals_db.type_lut);
set(gca,'XTickLabel',stateslabels_ttls);
title('Fraction of time on arousal state');
for ci = 1:length(animals_db.type_lut)
    X=Tstates(animals_db.type_list==ci &animals_db.arousal_cluster_list~=0, :);
    X=bsxfun(@rdivide, X, nansum(X,2));
    x{ci} = X;
end
for si = 1:length(stateslabels)
    [~,p(1)] = ttest2(x{1}(:,si), x{2}(:,si,:),'tail','right');
    [~,p(2)] = ttest2(x{1}(:,si), x{3}(:,si,:),'tail','right');
    disp(p)  
end


%% N states

% for ci = 1:length(animals_db.type_lut)
%     X=Nstates(animals_db.type_list==ci, :);
%     Mbytype(ci,:) = nanmean(X);
%     Sbytype(ci,:) = nanstd(X);
%     nbytype(ci,:) = sum(~isnan(X));
%
% end
% figure;b = barwitherr(Sbytype./sqrt(nbytype-1), Mbytype);
% for i=1:length(stateslabels_ttls)
%     b(i).FaceColor = CondColors(i,:);
% end
% set(gca,'XTickLabel', animals_db.type_lut);
% legend(stateslabels_ttls);
% title('Fraction of time on arousal state');
%
% mysave(gcf, fullfile(outputfiggolder, ['time_spent_fraction', num2str(length(stateslabels_ttls)) ]));


%% wheel histograms per state and behavior
binsN = linspace(0,3.5,100);

figure;
for ci = 1:length(animals_db.type_lut)
    animals = find(animals_db.type_list==ci &animals_db.arousal_cluster_list>0);
    subplot(1,length(animals_db.type_lut),ci);
    for state_i = 1:length(stateslabels)
        x = wheelvalsall.(stateslabels{state_i})(animals);
        
        set(gcf,'renderer','Painters');
        b=hist(x, binsN);
        bar(binsN,b/sum(b), 'facecolor',CondColors(state_i,:),'facealpha',.8,'edgecolor','none');
        title([ animals_db.type_lut{ci}]); xlabel('wheel speed');
        hold all
    end
    ylim([0 1]);
    legend(stateslabels_ttls)
end

mysave(gcf, fullfile(outputfiggolder, ['spont_wheel_hist_' num2str(length(stateslabels_ttls)) 'states']));

%% pupil histograms per state and behavior
figure;binsN = linspace(0,3.5,100);
for ci = 1:length(animals_db.type_lut)
    animals = find(animals_db.type_list==ci&animals_db.arousal_cluster_list>0);
    subplot(1,length(animals_db.type_lut),ci);
    for state_i = 1:length(stateslabels)
        x = pupvalsall.(stateslabels{state_i})(animals);
        
        set(gcf,'renderer','Painters');
        b=hist(x, binsN);
        bar(binsN,b/sum(b), 'facecolor',CondColors(state_i,:),'facealpha',.8,'edgecolor','none');
        title([ animals_db.type_lut{ci}]); xlabel('wheel speed');
        hold all
    end
    ylim([0 1]);
    legend(stateslabels_ttls)
end
mysave(gcf, fullfile(outputfiggolder, ['spont_pup_hist_' num2str(length(stateslabels_ttls)) 'states']));


for ci = 1:length(animals_db.type_lut)
    animals = find(animals_db.type_list==ci&animals_db.arousal_cluster_list>0);
    for state_i = 1:length(stateslabels)
        Mp(state_i,ci) = nanmean(pupvalsall.(stateslabels{state_i})(animals));
        Sp(state_i,ci) = nanstd(pupvalsall.(stateslabels{state_i})(animals));
        Np(state_i,ci) = sum(~isnan(pupvalsall.(stateslabels{state_i})(animals)));
        Mw(state_i,ci) = nanmean(wheelvalsall.(stateslabels{state_i})(animals));
        Sw(state_i,ci) = nanstd(wheelvalsall.(stateslabels{state_i})(animals));
        Nw(state_i,ci) = sum(~isnan(wheelvalsall.(stateslabels{state_i})(animals)));
    end
end
figure;subplot(2,1,1);b=barwitherr(Sp'./sqrt(Np'-1), Mp');
set(gca,'XTickLabel',animals_db.type_lut);
for si=1:length(stateslabels)
    b(si).FaceColor = CondColors(si,:);
end
ylabel('Pupil');legend(stateslabels_ttls);

subplot(2,1,2);b=barwitherr(Sw'./sqrt(Nw'-1), Mw');
set(gca,'XTickLabel',animals_db.type_lut);
for si=1:length(stateslabels)
    b(si).FaceColor = CondColors(si,:);
end
ylabel('Wheel');legend(stateslabels_ttls);
set(gca,'XTickLabel',animals_db.type_lut);
mysave(gcf, fullfile(outputfiggolder, ['mean_wheel_pupil', num2str(length(stateslabels_ttls)) ]));

end
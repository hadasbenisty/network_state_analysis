% Fig 1 - Behavior
function main_make_behavior_figs_crispr
animals_db = get_animals_meta_data_by_csv;
addpath(genpath('../functions'))
addpath(genpath('../utils'))
stateslabels = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
procfolder = 'X:\Hadas\Meso-imaging\crispr\meso_results\ProcessingDirectory\';
outputfiggolder = 'X:\Hadas\Meso-imaging\crispr\meso_results\figs\behavior';
mkNewDir(outputfiggolder);
mkNewDir(procfolder);
plot_time_spent(animals_db, stateslabels, procfolder, outputfiggolder)


end




function plot_time_spent(animals_db, stateslabels, procfolder, outputfiggolder)


if exist(fullfile(procfolder, 'time_spent_data.mat'), 'file')
    load(fullfile(procfolder, 'time_spent_data.mat'), 'wheelvalsall', 'pupvalsall');   
else
    n = length(animals_db.animal_list);
    for state_i = 1:length(stateslabels)
        wheelvalsall.(stateslabels{state_i})=cell(n,1);
        pupvalsall.(stateslabels{state_i})=cell(n,1);
        
    end
    
    Nstates = nan(n, length(stateslabels));
    for ai = 1:n
        load(fullfile(procfolder, animals_db.folder_list{ai}, 'arousal_state_ITI_segemts.mat'),...
            'segments_arousals','pupvals','wheelvals');
        for statei = 1:length(stateslabels)
            if isfield(segments_arousals, stateslabels{statei})
                Nstates(ai, statei) = sum(diff(segments_arousals.(stateslabels{statei}),1,2));
                wheelvalsall.(stateslabels{state_i}){ai} = cat(1,wheelvalsall.(stateslabels{state_i}){ai},...
                    wheelvals.(stateslabels{state_i}));
                pupvalsall.(stateslabels{state_i}){ai} = cat(1,pupvalsall.(stateslabels{state_i}){ai},...
                    pupvals.(stateslabels{state_i}));
            end
        end
    end
   save(fullfile(procfolder, 'time_spent_data.mat'),'Nstates', 'wheelvalsall', 'pupvalsall'); 
end
stateslabels_ttls = {'low pup q' 'high pup q' 'high pup l'};


%% time spent
figure;
for ci = 1:length(animals_db.type_lut)
    X=Nstates(animals_db.type_list==ci, :);
    X=bsxfun(@rdivide, X, nansum(X,2));
    Mbytype = nanmean(X);
    Sbytype = nanstd(X);
    nbytype = sum(~isnan(X));
    subplot(length(animals_db.type_lut),1,ci);
    plot_3_bars(Mbytype, Sbytype./sqrt(nbytype-1), stateslabels_ttls)
    ylim([0 1]);
    title(animals_db.type_lut{ci})
end
suptitle('Fraction of time on arousal state');

mysave(gcf, fullfile(outputfiggolder, 'time_spent_fraction'));


%% wheel histograms per state and behavior
binsN = linspace(0,60,32);

CondColors = get_3states_colors;
figure;l=1;    
for state_i = 1:length(stateslabels)
    for ci = 1:length(animals_db.type_lut)        
        animals = find(animals_db.type_list==ci);
        x=[];
        for k=1:length(animals)
            x = cat(1,x, wheelvalsall.(stateslabels{state_i}){animals(k)});
        end
        subplot(length(stateslabels),length(animals_db.type_lut),l);
        set(gcf,'renderer','Painters');
        histogram(x, binsN,'facecolor',CondColors(state_i,:),'facealpha',.8,'edgecolor','none');
        title([stateslabels_ttls{state_i} ' ' animals_db.type_lut{ci}]);xlim([0 60]); xlabel('wheel');
        l=l+1;
    end
end
mysave(gcf, fullfile(outputfiggolder, 'spont_wheel_hist_3states'));

%% pupil histograms per state and behavior
binsN = linspace(0,6e3,32);
figure;l=1;    
for state_i = 1:length(stateslabels)
    for ci = 1:length(animals_db.type_lut)        
        animals = find(animals_db.type_list==ci);
        x=[];
        for k=1:length(animals)
            x = cat(1,x, pupvalsall.(stateslabels{state_i}){animals(k)});
        end
        subplot(length(stateslabels),length(animals_db.type_lut),l);
        set(gcf,'renderer','Painters');
        histogram(x, binsN,'facecolor',CondColors(state_i,:),'facealpha',.8,'edgecolor','none');
        title([stateslabels_ttls{state_i} ' ' animals_db.type_lut{ci}]);xlim([0 60]); xlabel('wheel');
        l=l+1;
    end
end
mysave(gcf, fullfile(outputfiggolder, 'spont_pup_histogram_3states'));


end


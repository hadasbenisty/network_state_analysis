function main_make_network_figs_crispr
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
animals = get_animals_meta_data_by_csv;
stateslabels3 = { 'low_pupil', 'high_pupil', 'loc'};
stateslabels2 = { 'qui', 'loc'};
% cent_features = {  'eigenvector' 'degree' 'closeness' 'participation' , 'diffmap',  'betweenness' 'pagerank', 'second_eigval'};%};%'eigenvector'
cent_features = {  'eigenvector' 'degree'  'second_eigval'};%};%'eigenvector'
similarity_name = {'pearson_corr'   };%'corr',,  'fullcorr' 'cov''partial_corr'
doover=true;
%% Fig 4 - network
% plot_centrality_res_per_day(animals, outputfiggolder, stateslabels);
for sim_i = 1:length(similarity_name)
    
    outputfiggolder = ['X:\Hadas\Meso-imaging\CRISPR\Figures\network_centrality_' similarity_name{sim_i}];
    procfolder = ['x:\Hadas\Meso-imaging\CRISPR\analysis_results\network_centrality_' similarity_name{sim_i}];
    mkNewDir(outputfiggolder)
    
%     plot_correlation_matrics_permutation_test(procfolder, similarity_name{sim_i}, animals, outputfiggolder);
    
    plot_centrality_res(cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels3, doover);
    plot_centrality_across_groups_arousal(cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels3, doover);
    %     plot_ctrl_vs_mutants(cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels3, doover);
    
    plot_centrality_across_groups_arousal(cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels2, doover);
    %     plot_ctrl_vs_mutants(cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels2, doover);
    plot_centrality_res(cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels2, doover);
    %     plot_correlation_matrices_type_sex(stateslabels2, similarity_name{sim_i}, animals, outputfiggolder);
    
end
end
function plot_correlation_matrics_permutation_test(procfolder, simname, animals, outputfiggolder)
n=length(animals.folder_list);
validinds = animals.toinclude_list==find(strcmp(animals.toinclude_lut, 'Good'));
high_pupildata_sha = nan(23,23,1000,n);
locdata_sh = nan(23,23,1000,n);
highpupildata_shb = nan(23,23,1000,n);
lowpupildata_sh = nan(23,23,1000,n);
highfacedata_sh = nan(23,23,1000,n);
lowfacedata_sh = nan(23,23,1000,n);
sitdata = nan(23,23,n);
locdata = nan(23,23,n);
highpupildata = nan(23,23,n);
lowpupildata = nan(23,23,n);
highfacedata = nan(23,23,n);
lowfacedata = nan(23,23,n);

for ai = 1:length(animals.folder_list)
    if validinds(ai)
        if isfile(fullfile(procfolder,animals.folder_list{ai} ,'shuffled_corr_Allen_.mat'))
            
            shdata = load(fullfile(procfolder,animals.folder_list{ai} ,'shuffled_corr_Allen_.mat'));
            if isempty(shdata.corrmat_highpupilloc)||isempty(shdata.corrmat_highpupilloc.high_pupil)
                continue;
            end
            high_pupildata_sha(:,:,:,ai) = shdata.corrmat_highpupilloc.high_pupil;
            locdata_sh(:,:,:,ai) = shdata.corrmat_highpupilloc.loc;
            if isfield(shdata, 'corrmat_pupil')&&~isempty(shdata.corrmat_pupil)
                lowpupildata_sh(:,:,:,ai) = shdata.corrmat_pupil.low_pupil;
                highpupildata_shb(:,:,:,ai) = shdata.corrmat_pupil.high_pupil;
            end
           
        end
        if isfile(fullfile(procfolder,animals.folder_list{ai} ,'locAllen_Inf.mat'))
            load(fullfile(procfolder,animals.folder_list{ai} ,'locAllen_Inf.mat'),'W_corr');
            locdata(:,:,ai) = W_corr;
        end
        if isfile(fullfile(procfolder,animals.folder_list{ai} ,'sitAllen_Inf.mat'))
            load(fullfile(procfolder,animals.folder_list{ai} ,'sitAllen_Inf.mat'),'W_corr');
            sitdata(:,:,ai) = W_corr;
        end
        if isfile(fullfile(procfolder,animals.folder_list{ai} ,'high_pupilAllen_Inf.mat'))
            load(fullfile(procfolder,animals.folder_list{ai} ,'high_pupilAllen_Inf.mat'),'W_corr');
            highpupildata(:,:,ai) = W_corr;
        end
        if isfile(fullfile(procfolder,animals.folder_list{ai} ,'low_pupilAllen_Inf.mat'))
            load(fullfile(procfolder,animals.folder_list{ai} ,'low_pupilAllen_Inf.mat'),'W_corr');
            lowpupildata(:,:,ai) = W_corr;
        end
     
     
        
    end
end

for ci = 1:length(animals.type_lut)
    ii =  animals.type_list==ci&~squeeze(isnan(locdata(1,2,:)));
    nbytype(ci) = length(unique(animals.animal_list(ii)));
end
for ci = 1:length(animals.type_lut)
     % 2->3
    [hpupil23(:,:,ci),pvalpupil23(:,:,ci)]=corr_perm_test(locdata, highpupildata, locdata_sh, high_pupildata_sha, animals.type_list==ci);
    % 1->2
    [hpupil12(:,:,ci),pvalpupil12(:,:,ci)]=corr_perm_test(highpupildata, lowpupildata, highpupildata_shb, lowpupildata_sh, animals.type_list==ci);
   
end
% 1->2
plot_corr_by_state(animals, highpupildata, lowpupildata, hpupil12, 'High Pupil', 'Low Pupil');
mysave(gcf, fullfile(outputfiggolder, 'corr_pupil_1_2'));
% 2->3
plot_corr_by_state(animals, locdata, highpupildata, hpupil23, 'Loc', 'High Pupil');
mysave(gcf, fullfile(outputfiggolder, 'corr_pupil_2_3'));



end
function [trial_states, spon_states] = load_similarity_matrices(signame, outputfolder, animals, statenames, isweigtedstr)

for i=1:length(animals)
    animal=animals{i};
    for state_i = 1:length(statenames)
        % trial correct
        load(fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame]), ...
            'W_corr_cor');
        load(fullfile(outputfolder,[animal '_',statenames{state_i} ,'trials_incorrect_' signame]), ...
            'W_corr_inc');
        load(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signame]), ...
            'W_corr');
        switch signame
            case 'Allen'
                
                
                spon_states(:,:,state_i, i) = process_sim(W_corr); %#ok<AGROW>
                trial_states.correct(:,:,state_i, i) = process_sim(W_corr_cor);
                trial_states.incorrect(:,:,state_i, i) = process_sim(W_corr_inc);
            case 'LSSC'
                spon_states.(statenames{state_i}){i}  = process_sim(W_corr);
                trial_states.(statenames{state_i}).correct{i} = process_sim(W_corr_cor);
                trial_states.(statenames{state_i}).incorrect{i} = process_sim(W_corr_inc);
        end
        
    end
end


end

function [trial_states, spon_states, spont_heatmap, trial_heatmap]  = load_centrality_results_per_days_allen(cent_features, signame, outputfolder, animals, statenames, isweigtedstr)
for state_i = 1:length(statenames)
    for cent_i = 1:length(cent_features)
        spon_states.(statenames{state_i}).(cent_features{cent_i}) = nan(23, 40, 6);
        trial_states.(statenames{state_i}).(cent_features{cent_i}).correct = nan(23, 40, 6);
        trial_states.(statenames{state_i}).(cent_features{cent_i}).incorrect =nan(23, 40, 6);
        spont_heatmap.(statenames{state_i}).(cent_features{cent_i}) = zeros(256, 256, length(animals));
        trial_heatmap.(statenames{state_i}).(cent_features{cent_i}).correct = zeros(256, 256, length(animals));
        trial_heatmap.(statenames{state_i}).(cent_features{cent_i}).incorrect = zeros(256, 256, length(animals));
        
    end
end


for i=1:length(animals)
    animal=animals{i};
    [~,days_list] = animaltodays(animal);
    for day_i = 1:length(days_list)
        
        for state_i = 1:length(statenames)
            % trial correct
            filename = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_day' num2str(days_list(day_i)) '.mat']);
            if exist(filename, 'file')
                load(filename, ['cent_corr_' isweigtedstr ]);
                centvals = eval(['cent_corr_' isweigtedstr ]);
                for cent_i = 1:length(cent_features)
                    switch signame
                        case 'Allen'
                            trial_states.(statenames{state_i}).(cent_features{cent_i}).correct(:, days_list(day_i), i) = centvals.(cent_features{cent_i});
                    end
                end
            end
            % trial incorrect
            filename = fullfile(outputfolder,[animal '_',statenames{state_i} ,'trials_incorrect_' signame '_day' num2str(days_list(day_i)) '.mat']);
            if exist(filename, 'file')
                load(filename, ['cent_corr_' isweigtedstr ]);
                centvals = eval(['cent_corr_' isweigtedstr ]);
                for cent_i = 1:length(cent_features)
                    switch signame
                        case 'Allen'
                            trial_states.(statenames{state_i}).(cent_features{cent_i}).incorrect(:, days_list(day_i), i) = centvals.(cent_features{cent_i});
                    end
                end
                
            end
            % spont
            filename = fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signame '_day' num2str(days_list(day_i)) '.mat']);
            if exist(filename, 'file')
                load(filename, ['cent_corr_' isweigtedstr ]);
                centvals = eval(['cent_corr_' isweigtedstr ]);
                for cent_i = 1:length(cent_features)
                    switch signame
                        case 'Allen'
                            spon_states.(statenames{state_i}).(cent_features{cent_i})(:, days_list(day_i), i) = centvals.(cent_features{cent_i});
                    end
                    
                end
            end
        end
    end
end
end

function [spon_states, spont_heatmap]  = load_centrality_results(cent_features, signame, outputfolder, animals, statenames, isweigtedstr, th)

spon_states = nan(23, length(animals), length(statenames), length(cent_features));
spont_heatmap = nan(256,256, length(animals),length(statenames), length(cent_features));






for i=1:length(animals)
    animal=animals{i};
    for state_i = 1:length(statenames)
        
        if ~exist(fullfile(outputfolder,animal,[statenames{state_i} ,signame '_' num2str(th) '.mat']), 'file')
            continue;
        end
        
        % spont
        load(fullfile(outputfolder,animal,[statenames{state_i} ,signame '_' num2str(th)]), ...
            ['cent_corr_' isweigtedstr ] );
        centvals = eval(['cent_corr_' isweigtedstr ]);
        
        
        xx=[];
        for cent_i = 1:length(cent_features)
            xx(:, cent_i) = centvals.(cent_features{cent_i});
            
        end
        Pvec = scores_to_heatmap(xx, 0,signame, animal);
        spon_states(:, i, state_i, :) = heatmap2score_allen(Pvec);
        for cent_i = 1:length(cent_features)
            spont_heatmap(:,:,i, state_i,cent_i)=Pvec(:,:,cent_i);
        end
        
        
    end
    
end
end

function plot_correlation_matrices_type_sex(statenames, simname, animals, outputfiggolder)

outputfolder=['X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr\network_centrality_' simname];
files = dir('X:\Hadas\Meso-imaging\lan\xxpsych\spt\xx_12_grid4.mat');
load(fullfile(files(1).folder, files(1).name), 'par_inds');

signames = {'Allen' };%

parcels_names = get_allen_meta_parcels;


l=1;labels_type=[];labels_state=[];
for state_i = 1:length(statenames)
    for ai = 1:length(animals.folder_list)
        if animals.toinclude_list(ai)==find(strcmp(animals.toinclude_lut, 'Good'))
            filename=fullfile(outputfolder,animals.folder_list{ai} ,[statenames{state_i} ,signames{1} '_Inf.mat']);
            if ~isfile(filename)
                continue;
            end
            r= load(filename,'W_corr','cent_corr_weighted');
            
            Wvec(:,:,l) = r.W_corr+eye(23);
            deg(:,l) = r.cent_corr_weighted.degree;
            eigv(:,l) = r.cent_corr_weighted.eigenvector;
            labels_type(l) = animals.type_list(ai);
            labels_state(l) = state_i;
            animal_label(l) = ai;
            cohort_label(l) = animals.cohort_list(ai);
            sex_label(l) = animals.sex_list(ai);
            
            l=l+1;
        end
    end
end


%% by sex and type
figure;
for sex_i = 2:3
    M=[];S=[];
    for type_i = 1:length(animals.type_lut)
        
        M(:,type_i) = mean(eigv(:,labels_type==type_i),2);
        S(:,type_i) = std(eigv(:,labels_type==type_i),[],2)./sqrt(sum(labels_type==type_i)-1);
    end
    subplot(2,1,sex_i-1);barwitherr(S,M);set(gca,'XTick',1:length(parcels_names));
    set(gca,'XTickLabel',parcels_names);title(animals.sex_lut{sex_i});legend(animals.type_lut)
end
mysave(gcf,fullfile(outputfiggolder,'eigenvector','by_sex'));

for co=1:3
    figure;i=1;
    for sex_i = 2:3
        for type_i = 1:length(animals.type_lut)
            subplot(2,length(animals.type_lut),i);i=i+1;
            imagesc(mean(Wvec(:,:,cohort_label==co&labels_type==type_i&sex_label==sex_i),3));
            title([animals.type_lut{type_i} ' ' animals.sex_lut{sex_i}]);
            set(gca,'XTick',[]);
            set(gca,'YTick',1:2:length(parcels_names));
            set(gca,'YTickLabel',parcels_names(1:2:end));
        end
    end
    suptitle(animals.cohort_lut{co});
    mysave(gcf,fullfile(outputfiggolder, ['corr_by_sex_type_co' num2str(co)]));
end



end

function plot_centrality_across_groups_arousal(cent_features, simname, animals, outputfiggolder, statenames, doover)

outputfolder=['X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr\network_centrality_' simname];
files = dir(['X:\Hadas\Meso-imaging\lan\xxpsych\spt\xx_12_grid4.mat']);
load(fullfile(files(1).folder, files(1).name), 'par_inds');

signals_names = {'Allen' };%
isweigtedstr = { 'weighted'};%'notweighted'
for l = 1:length(cent_features)
    mkNewDir(fullfile(outputfiggolder,cent_features{l}));
end
[parcels_names, parcels_region_labels, final_index, region_lut, grid_map_final_index, labelsbyallen] = get_allen_meta_parcels;
legstr = statenames;
for k=1:length(legstr)
    legstr{k}(legstr{k}=='_') = ' ';
end

for sig_i = 1:length(signals_names)
    switch signals_names{sig_i}
        case 'Allen'
            visualinds = 1;
            somatoinds = 14;
            visualinds = find(parcels_region_labels==1);
            somatoinds = find(parcels_region_labels==6);
            thT=[ Inf ];
        case 'Grid4'
            %             [parcels_names, parcels_region_labels, final_index_grid, region_lut, grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(par_inds, 4);
            %             visualinds = find(parcels_region_labels==1);
            %             somatoinds = find(parcels_region_labels==6);
            thT=[ Inf ];
    end
    arousaltypes = unique(animals.type_list);
    Nstates = length(statenames);
    for th=thT
        for ti = 1:length(arousaltypes)
            curtype = animals.type_lut{ti};
            animalsinds = find(animals.type_list==ti & animals.toinclude_list == find(strcmp(animals.toinclude_lut, 'Good'))&animals.arousal_cluster_list>0);
            for isweigted = 1:length(isweigtedstr)
                
                sumfile = fullfile(outputfolder, [num2str(Nstates) 'states_' curtype '_summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'th' num2str(th) '.mat']);
                if  ~doover&&exist(sumfile, 'file')
                    load(sumfile);
                else
                    [spon_states, spont_heatmap] = load_centrality_results(cent_features, signals_names{sig_i}, outputfolder, animals.folder_list(animalsinds), statenames, isweigtedstr{isweigted}, th);
                    save(sumfile, 'spon_states',...
                        'spont_heatmap');
                end
                cent_by_allen_mean(:,:,:,ti) = squeeze(nanmean(spon_states,2));
                cent_by_allen_std(:,:,:,ti) = squeeze(nanstd(spon_states,[],2));
                cent_by_allen_N(:,:,ti) = squeeze(sum(~isnan(spon_states(1,:,:,:,:)),2));
                masks_cent_mean_spont(:,:, :, :,ti) = squeeze(nanmean(spont_heatmap,3));
                %% corr matrices
                [~, W_spont(:,:,:,ti)]=plot_corr_matrices_allen_crispr(outputfolder, animals.folder_list(animalsinds),statenames, th, parcels_names, parcels_region_labels);
                close all;
            end
        end
        if ~all(isnan(W_spont(:)))
            
            %% type & arousal state
            figure;l=1;
            for si=1:length(statenames)
                for ti=1:length(arousaltypes)
                    subplot(length(statenames)+1,length(arousaltypes),l);
                    plot_corr_mat(W_spont(:,:,si,ti),parcels_names,[0 1]);
                    title([ animals.arousal_cluster_lut{2+ti} ' ' statenames{si} ]);l=l+1;
                    colorbar;
                end
            end
            colormap jet;
            for ti=1:length(arousaltypes)
                subplot(length(statenames)+1,length(arousaltypes),l);
                plot_corr_mat(W_spont(:,:,end,ti)-W_spont(:,:,1,ti),parcels_names,[-.5  .5]/5);l=l+1;
                title([ animals.arousal_cluster_lut{2+ti} ' ' statenames{end} '-' statenames{1} ]);
                colorbar; %
            end
            colormap(redblue);
            mysave(gcf, fullfile(outputfiggolder,[num2str(Nstates) 'states_All_types_spont_', '_W_' signals_names{sig_i} '_th' num2str(th)]));
        end
        
        
        figure;l=1;
        for si=1:length(statenames)
            for ti=2:length(arousaltypes)
                subplot(length(statenames),length(arousaltypes)-1,l);
                plot_corr_mat(W_spont(:,:,si,1)-W_spont(:,:,si,ti),parcels_names,[-0.2 .2]);
                l=l+1;
                title([ 'Ctrl minus ' animals.type_lut{ti} ' ' statenames{si} ]);
                colorbar;colormap(redblue);
            end
        end
        mysave(gcf, fullfile(outputfiggolder,[num2str(Nstates) 'statesCorrMatrices_Ctrlminus_types_spont_', '_W_' signals_names{sig_i} '_th' num2str(th)]));
        
        
        ni = find(strcmp(cent_features, 'second_eigval'));
        M = squeeze(cent_by_allen_mean(1,:,ni,:));
        S = squeeze(cent_by_allen_std(1,:,ni,:));
        N = squeeze(cent_by_allen_N(:,ni,:));
        figure;
        barwitherr(S./sqrt(N-1),M);
        legend(animals.type_lut);
        set(gca,'XTickLabel',(legstr));
        title('Second Eigenval');
        mysave(gcf, fullfile(outputfiggolder,'second_eigval',[num2str(Nstates) 'states_All_types_spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
        
        
        %         [spon_states, spont_heatmap] = load_centrality_results(cent_features, signals_names{sig_i}, outputfolder, animals.folder_list, statenames, isweigtedstr{isweigted}, th);
        
        for ni = find(~strcmp(cent_features, 'second_eigval'))
            %             x=[];M=[];N=[];S=[];
            %             for ti=1:length(statenames)
            %                 x{ti} = spon_states(:,animals.arousal_cluster_list>0&animals.type_list==ti,end,ni)-...
            %                     spon_states(:,animals.arousal_cluster_list>0&animals.type_list==ti,1,ni);
            %                 ii = ~isnan(x{ti}(1,:));
            %                 x{ti}=x{ti}(:,ii);
            %                 M(:,ti) = nanmean(x{ti},2);
            %                 for r=1:100
            %                     R = randn([size(x{ti},2),1])>0;
            %                     xs{ti} = -x{ti}(:,R==1);
            %                     Ms(:,ti,r) = nanmean(xs{ti},2);
            %                 end
            %             end
            %             Ms=nanmean(Ms,3); figure;
            %             for ti=1:length(statenames)
            %            subplot(3,1,ti);bar([M(:,ti) Ms(:,ti)]);
            %            legend('Data','Shuffle');
            %             end
            %             set(gca,'XTick',1:23);set(gca,'XTickLabel',parcels_names)
            %
            %
            M = permute(squeeze(cent_by_allen_mean(:,:,ni,:)),[1 3 2]);
            S = permute(squeeze(cent_by_allen_std(:,:,ni,:)),[1 3 2]);
            N = repmat(permute(squeeze(cent_by_allen_N(:,ni,:)),[3 2 1]), 23, 1,1);
            
            legend(animals.type_lut);set(gca,'XTick',1:23);set(gca,'XTickLabel',parcels_names)
            for si=1:length(statenames)
                h(si)=subplot(length(statenames),1, si);
                barwitherr(S(:,:,si)./sqrt(N(:,:,si)-1),M(:,:,si));
                legend(animals.type_lut);
                set(gca,'XTick',1:length(parcels_names))
                set(gca,'XTickLabel',(parcels_names));
                title(legstr{si});
            end
            linkaxes(h,'y');
            suptitle(cent_features{ni});
            mysave(gcf, fullfile(outputfiggolder,cent_features{ni}, [num2str(Nstates) 'states_All_types_spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
            
            L = quantile(reshape(masks_cent_mean_spont(:,:,:,ni,:),1,[]),[.1 .9]);
            figure;l=1;
            for si=1:length(statenames)
                for ti=1:length(arousaltypes)
                    subplot(length(statenames)+1,length(arousaltypes),l);
                    plot_vals_heatmap(masks_cent_mean_spont(:,:,si,ni,ti),...
                        '',[],  L(1), L(2), 1,colormap(redblue));l=l+1;
                    title([ animals.type_lut{ti} ' ' statenames{si} ]);
                end
            end
            
            L = quantile(reshape(masks_cent_mean_spont(:,:,end,ni,:)-masks_cent_mean_spont(:,:,1,ni,:),1,[]),[.1 .9]);
            for ti=1:length(arousaltypes)
                subplot(length(statenames)+1,length(arousaltypes),l);
                plot_vals_heatmap(masks_cent_mean_spont(:,:,end,ni,ti)-masks_cent_mean_spont(:,:,1,ni,ti),...
                    '',[],  -max(abs(L)), max(abs(L)), 1,colormap(redblue));l=l+1;
                title([ animals.type_lut{ti} ' ' statenames{end} '-' statenames{1} ]);
            end
            suptitle(cent_features{ni});
            mysave(gcf, fullfile(outputfiggolder,cent_features{ni}, [num2str(Nstates) 'states_All_types_spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_heatmaps_' signals_names{sig_i} '_th' num2str(th)]));
            
        end
        
    end
end

end
function plot_centrality_across_clusters_arousal(cent_features, simname, animals, outputfiggolder, statenames, doover)

% outputfolder=['X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr\network_centrality_' simname];
files = dir(['X:\Hadas\Meso-imaging\lan\xxpsych\spt\xx_12_grid4.mat']);
load(fullfile(files(1).folder, files(1).name), 'par_inds');

signals_names = {'Allen' };%
isweigtedstr = { 'weighted'};%'notweighted'
for l = 1:length(cent_features)
    mkNewDir(fullfile(outputfiggolder,cent_features{l}));
end
[parcels_names, parcels_region_labels, final_index, region_lut, grid_map_final_index, labelsbyallen] = get_allen_meta_parcels;
legstr = statenames;
for k=1:length(legstr)
    legstr{k}(legstr{k}=='_') = ' ';
end

for sig_i = 1:length(signals_names)
    switch signals_names{sig_i}
        case 'Allen'
            visualinds = 1;
            somatoinds = 14;
            visualinds = find(parcels_region_labels==1);
            somatoinds = find(parcels_region_labels==6);
            thT=[ Inf ];
        case 'Grid4'
            %             [parcels_names, parcels_region_labels, final_index_grid, region_lut, grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(par_inds, 4);
            %             visualinds = find(parcels_region_labels==1);
            %             somatoinds = find(parcels_region_labels==6);
            thT=[ Inf ];
    end
    arousaltypes = 1:3;
    Nstates = length(statenames);
    for th=thT
        for ti = 1:length(arousaltypes)
            curtype = animals.type_lut{ti};
            animalsinds = find(animals.toinclude_list == find(strcmp(animals.toinclude_lut, 'Good')));
            for isweigted = 1:length(isweigtedstr)
                
                sumfile = fullfile(outputfolder, [num2str(Nstates) 'states_' curtype '_summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'th' num2str(th) '.mat']);
                if  ~doover&&exist(sumfile, 'file')
                    load(sumfile);
                else
                    [spon_states, spont_heatmap] = load_centrality_results(cent_features, signals_names{sig_i}, outputfolder, animals.folder_list(animalsinds), statenames, isweigtedstr{isweigted}, th);
                    save(sumfile, 'spon_states',...
                        'spont_heatmap');
                end
                cent_by_allen_mean(:,:,:,ti) = squeeze(nanmean(spon_states,2));
                cent_by_allen_std(:,:,:,ti) = squeeze(nanstd(spon_states,[],2));
                cent_by_allen_N(:,:,ti) = squeeze(sum(~isnan(spon_states(1,:,:,:,:)),2));
                masks_cent_mean_spont(:,:, :, :,ti) = squeeze(nanmean(spont_heatmap,3));
                %% corr matrices
                [~, W_spont(:,:,:,ti)]=plot_corr_matrices_allen_crispr(outputfolder, animals.folder_list(animalsinds),statenames, th, parcels_names, parcels_region_labels);
                close all;
            end
        end
        if ~all(isnan(W_spont(:)))
            
            
            %% type & arousal state
            figure;l=1;
            for si=1:length(statenames)
                for ti=1:length(arousaltypes)
                    subplot(length(statenames)+1,length(arousaltypes),l);
                    plot_corr_mat(W_spont(:,:,si,ti),parcels_names,[0 1]);
                    title([ animals.arousal_cluster_lut{2+ti} ' ' statenames{si} ]);l=l+1;
                    colorbar;
                end
            end
            colormap jet;
            for ti=1:length(arousaltypes)
                subplot(length(statenames)+1,length(arousaltypes),l);
                plot_corr_mat(W_spont(:,:,end,ti)-W_spont(:,:,1,ti),parcels_names,[-.5  .5]/5);l=l+1;
                title([ animals.arousal_cluster_lut{2+ti} ' ' statenames{end} '-' statenames{1} ]);
                colorbar; %
            end
            colormap(redblue);
            mysave(gcf, fullfile(outputfiggolder,[num2str(Nstates) 'states_All_arousalclusters_spont_', '_W_' signals_names{sig_i} '_th' num2str(th)]));
        end
        
        
        figure;l=1;
        for si=1:length(statenames)
            for ti=2:length(arousaltypes)
                subplot(length(statenames),length(arousaltypes)-1,l);
                plot_corr_mat(W_spont(:,:,si,1)-W_spont(:,:,si,ti),parcels_names,[-0.2 .2]);
                l=l+1;
                title([ 'Ctrl minus ' animals.arousal_cluster_lut{2+ti} ' ' statenames{si} ]);
                colorbar;colormap(redblue);
            end
        end
        mysave(gcf, fullfile(outputfiggolder,[num2str(Nstates) 'statesCorrMatrices_Ctrlminus__arousalclusters_spont_', '_W_' signals_names{sig_i} '_th' num2str(th)]));
        
        
        ni = find(strcmp(cent_features, 'second_eigval'));
        M = squeeze(cent_by_allen_mean(1,:,ni,:));
        S = squeeze(cent_by_allen_std(1,:,ni,:));
        N = squeeze(cent_by_allen_N(:,ni,:));
        figure;
        barwitherr(S./sqrt(N-1),M);
        legend(animals.type_lut);
        set(gca,'XTickLabel',(legstr));
        title('Second Eigenval');
        mysave(gcf, fullfile(outputfiggolder,'second_eigval',[num2str(Nstates) 'states_All_types_spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
        
        
        %         [spon_states, spont_heatmap] = load_centrality_results(cent_features, signals_names{sig_i}, outputfolder, animals.folder_list, statenames, isweigtedstr{isweigted}, th);
        
        for ni = find(~strcmp(cent_features, 'second_eigval'))
            %             x=[];M=[];N=[];S=[];
            %             for ti=1:length(statenames)
            %                 x{ti} = spon_states(:,animals.arousal_cluster_list>0&animals.type_list==ti,end,ni)-...
            %                     spon_states(:,animals.arousal_cluster_list>0&animals.type_list==ti,1,ni);
            %                 ii = ~isnan(x{ti}(1,:));
            %                 x{ti}=x{ti}(:,ii);
            %                 M(:,ti) = nanmean(x{ti},2);
            %                 for r=1:100
            %                     R = randn([size(x{ti},2),1])>0;
            %                     xs{ti} = -x{ti}(:,R==1);
            %                     Ms(:,ti,r) = nanmean(xs{ti},2);
            %                 end
            %             end
            %             Ms=nanmean(Ms,3); figure;
            %             for ti=1:length(statenames)
            %            subplot(3,1,ti);bar([M(:,ti) Ms(:,ti)]);
            %            legend('Data','Shuffle');
            %             end
            %             set(gca,'XTick',1:23);set(gca,'XTickLabel',parcels_names)
            %
            %
            M = permute(squeeze(cent_by_allen_mean(:,:,ni,:)),[1 3 2]);
            S = permute(squeeze(cent_by_allen_std(:,:,ni,:)),[1 3 2]);
            N = repmat(permute(squeeze(cent_by_allen_N(:,ni,:)),[3 2 1]), 23, 1,1);
            
            %             legend(animals.type_lut);set(gca,'XTick',1:23);set(gca,'XTickLabel',parcels_names)
            figure;
            for si=1:length(statenames)
                h(si)=subplot(length(statenames),1, si);
                barwitherr(S(:,:,si)./sqrt(N(:,:,si)-1),M(:,:,si));
                legend(animals.type_lut);
                set(gca,'XTick',1:length(parcels_names))
                set(gca,'XTickLabel',(parcels_names));
                title(legstr{si});
            end
            linkaxes(h,'y');
            suptitle(cent_features{ni});
            mysave(gcf, fullfile(outputfiggolder,cent_features{ni}, [num2str(Nstates) 'states_All_arousalclusters_spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
            
            L = quantile(reshape(masks_cent_mean_spont(:,:,:,ni,:),1,[]),[.1 .9]);
            figure;l=1;
            for si=1:length(statenames)
                for ti=1:length(arousaltypes)
                    subplot(length(statenames)+1,length(arousaltypes),l);
                    plot_vals_heatmap(masks_cent_mean_spont(:,:,si,ni,ti),...
                        '',[],  L(1), L(2), 1,colormap(redblue));l=l+1;
                    title([animals.arousal_cluster_lut{2+ti} ' ' statenames{si} ]);
                end
            end
            
            L = quantile(reshape(masks_cent_mean_spont(:,:,end,ni,:)-masks_cent_mean_spont(:,:,1,ni,:),1,[]),[.1 .9]);
            for ti=1:length(arousaltypes)
                subplot(length(statenames)+1,length(arousaltypes),l);
                plot_vals_heatmap(masks_cent_mean_spont(:,:,end,ni,ti)-masks_cent_mean_spont(:,:,1,ni,ti),...
                    '',[],  -max(abs(L)), max(abs(L)), 1,colormap(redblue));l=l+1;
                title([ animals.arousal_cluster_lut{2+ti} ' ' statenames{end} '-' statenames{1} ]);
            end
            suptitle(cent_features{ni});
            mysave(gcf, fullfile(outputfiggolder,cent_features{ni}, [num2str(Nstates) 'states_All_arousalclusters_spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_heatmaps_' signals_names{sig_i} '_th' num2str(th)]));
            
        end
        
    end
end

end













function plot_ctrl_vs_mutants(cent_features, simname, animals, outputfiggolder, statenames, doover)

% outputfolder=['X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr\network_centrality_' simname];
files = dir(['X:\Hadas\Meso-imaging\lan\xxpsych\spt\xx_12_grid4.mat']);
load(fullfile(files(1).folder, files(1).name), 'par_inds');

signals_names = {  'Allen'};%
thvals = [inf inf];
isweigtedstr = { 'weighted'};%'notweighted'
for l = 1:length(cent_features)
    mkNewDir(fullfile(outputfiggolder,cent_features{l}));
end
[parcels_names, parcels_region_labels, final_index, region_lut, grid_map_final_index, labelsbyallen] = get_allen_meta_parcels;

for sig_i = 1:length(signals_names)
    switch signals_names{sig_i}
        case 'Allen'
            visualinds = 1;
            somatoinds = 14;
            visualinds = find(parcels_region_labels==1);
            somatoinds = find(parcels_region_labels==6);
            thT=[ Inf ];%7
        case 'Grid4'
            %             [parcels_names, parcels_region_labels, final_index_grid, region_lut, grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(par_inds, 4);
            %             visualinds = find(parcels_region_labels==1);
            %             somatoinds = find(parcels_region_labels==6);
            thT=[ Inf ];%150
    end
    arousaltypes = unique(animals.type_list);
    Nstates = length(statenames);
    for th=thT
        for ti = 1:length(arousaltypes)
            curtype = animals.type_lut{ti};
            animalsinds = find(animals.type_list==ti & animals.toinclude_list==3&animals.arousal_cluster_list>0);
            isweigted = 1;
            
            sumfile = fullfile(outputfolder, [num2str(Nstates) 'states_' curtype '_summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'th' num2str(th) '.mat']);
            if  ~doover&&exist(sumfile, 'file')
                load(sumfile);
            else
                [spon_states, spont_heatmap] = load_centrality_results(cent_features, signals_names{sig_i}, outputfolder, animals.folder_list(animalsinds), statenames, isweigtedstr{isweigted}, th);
                save(sumfile, 'spon_states',...
                    'spont_heatmap');
            end
            
            
            maskscentmeanspont = squeeze(nanmean(spont_heatmap,3));
            ni = 1;
            for state_i=1:2
                Map(:,:,ti, state_i) =   nanmean(maskscentmeanspont(:,:,state_i,ni),4);
            end
        end
        
        
        for ti = 2:3
            for state_i=1:2
                figure;
                diffvals = Map(:,:,1, state_i)-Map(:,:,ti, state_i);
                
                L = quantile(diffvals(:),[.1 .9]);
                % spont
                %         subplot(2,3,3+ti);
                plot_vals_heatmap(diffvals, 'Centrality',...
                    [],  -abs(L(1)), abs(L(1)), 1,colormap(redblue))
                title(['Control-' animals.type_lut{ti} ' ' statenames{state_i}]);
                mysave(gcf,['ctrl_minus_' animals.type_lut{ti} '_' statenames{state_i} '_eigenvector']);
            end
            
        end
    end
end
end

function plot_centrality_res(cent_features, simname, animals, outputfiggolder, statenames, doover)

 outputfolder=['X:\Hadas\Meso-imaging\CRISPR\analysis_results\network_centrality_' simname];
files = dir(['X:\Hadas\Meso-imaging\lan\xxpsych\spt\xx_12_grid4.mat']);
load(fullfile(files(1).folder, files(1).name), 'par_inds');

signals_names = {'Allen'  };%'Grid4'
isweigtedstr = { 'weighted'};%'notweighted'
for l = 1:length(cent_features)
    mkNewDir(fullfile(outputfiggolder,cent_features{l}));
end
[parcels_names, parcels_region_labels, final_index, region_lut, grid_map_final_index, labelsbyallen] = get_allen_meta_parcels;

for sig_i = 1:length(signals_names)
    switch signals_names{sig_i}
        case 'Allen'
            visualinds = 1;
            somatoinds = 14;
            visualinds = find(parcels_region_labels==1);
            somatoinds = find(parcels_region_labels==6);
            thT=[ Inf ];%7
        case 'Grid4'
            %             [parcels_names, parcels_region_labels, final_index_grid, region_lut, grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(par_inds, 4);
            %             visualinds = find(parcels_region_labels==1);
            %             somatoinds = find(parcels_region_labels==6);
            thT=[ Inf ];%150
    end
    arousaltypes = unique(animals.type_list);
    Nstates = length(statenames);
    for ti = 1:length(arousaltypes)
        curtype = animals.type_lut{ti};
        animalsinds = find(animals.type_list==ti & animals.toinclude_list==3);
        for isweigted = 1:length(isweigtedstr)
            for th=thT
                sumfile = fullfile(outputfolder, [num2str(Nstates) 'states_' curtype '_summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'th' num2str(th) '.mat']);
                if  ~doover&&exist(sumfile, 'file')
                    load(sumfile);
                else
                    [spon_states, spont_heatmap] = load_centrality_results(cent_features, signals_names{sig_i}, outputfolder, animals.folder_list(animalsinds), statenames, isweigtedstr{isweigted}, th);
                    save(sumfile, 'spon_states',...
                        'spont_heatmap');
                end
                
                %% corr matrices
                %                 N=plot_corr_matrices_allen_crispr(outputfolder, animals.folder_list(animalsinds),statenames, th, parcels_names, parcels_region_labels);
                %                 simnames=simname;
                %                 simnames(simname=='_')=' ';
                %                 suptitle(['N=' num2str(N) ' ' simnames ' by states k=' num2str(th)]);
                %                 set(gcf,'Position',[680         104        1107         874]);
                %                 mysave(gcf, fullfile(outputfiggolder,[num2str(Nstates) 'states_' 'W_'   signals_names{sig_i} '_th' num2str(th) '_' curtype]));
                
                
                diffmapind = find(strcmp(cent_features, 'diffmap'));
                if ~isempty(diffmapind)
                    for j=1:length(statenames)
                        spon_states(:, :, j, diffmapind) = sign(diffmap2clusters(spon_states(:, :, j, diffmapind)));
                    end
                end
                
                legstr = statenames;
                for k=1:length(legstr)
                    legstr{k}(legstr{k}=='_') = ' ';
                end
                parcels_names = get_allen_meta_parcels;
                
                for ni = 1:2%find(~strcmp(cent_features, 'second_eigval'))
                    %                     for ii=1:length(animalsinds)
                    %                         if all(all(isnan(squeeze(spon_states(:,ii,:,ni)))))
                    %                             continue;
                    %                         end
                    %                         plot_bars_by_conditions(spon_states(:,ii,:,ni), ...
                    %                             cent_features{ni}, parcels_names, legstr, curtype);
                    %                         str = animals.folder_list{animalsinds(ii)};str(strfind(str, '/')) = '_';str(strfind(str, '\')) = '_';
                    %                         mysave(gcf, fullfile(outputfiggolder,cent_features{ni}, [str num2str(Nstates) 'states_' curtype '_spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
                    %                     end
                    plot_bars_by_conditions(spon_states(:,:,:,ni), ...
                        cent_features{ni}, parcels_names, legstr, curtype);
                    mysave(gcf, fullfile(outputfiggolder,cent_features{ni}, [num2str(Nstates) 'states_' curtype '_spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
                    
                end
                
                
                %% heatmaps
                
                maskscentmeanspont = squeeze(nanmean(spont_heatmap,3));
                for ni = 1:2%find(~strcmp(cent_features, 'second_eigval'))
                    
                    
                    figure;
                    v=reshape(maskscentmeanspont(:,:,:,ni),1,[]);
                    L = quantile(v,[.1 .9]);
                    if L(1)==L(2)
                        continue;
                    end
                    for state_i = 1:length(statenames)
                        % spont
                        subplot(2,2,state_i);
                        
                        plot_vals_heatmap(maskscentmeanspont(:,:,state_i,ni), 'Centrality',...
                            [],  L(1), L(2), 1,colormap(redblue))
                        title([statenames{state_i} ' Spont']);
                    end
                    diffmask = maskscentmeanspont(:,:,end,ni)-maskscentmeanspont(:,:,1,ni);
                    L = quantile(diffmask(:),[.1 .9]);
                    subplot(2,2,4);
                    plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
                        [],  -abs(L(1)), abs(L(1)), 1,colormap(redblue))
                    title([legstr{end} '-' legstr{1} ' spont']);
                    
                    
                    
                    suptitle(cent_features{ni});
                    mysave(gcf, fullfile(outputfiggolder,cent_features{ni},[num2str(Nstates) 'states_' curtype '_' cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i} '_th' num2str(th)]));
                    
                    
                    
                    
                    
                    
                    
                end
                %                 secondegvali = find(strcmp(cent_features, 'second_eigval'));
                %                 if ~isempty(secondegvali)
                %                     M = squeeze(nanmean(spon_states(1,:,:,secondegvali),2));
                %                     N = squeeze(sum(~isnan(spon_states(1,:,:,secondegvali)),2));
                %                     S = squeeze(nanstd(spon_states(1,:,:,secondegvali),[],2))./sqrt(N-1);
                %
                %                     figure;plot_2_bars(M,S,statenames);
                %
                %                     mysave(gcf, fullfile(outputfiggolder,[num2str(Nstates) 'states_' curtype '_second_eigval_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
                %                 end
                close all;
            end
            
        end
    end
end
end
% spont
%     diffmask = scores_to_heatmap_allen(mean(spon_states_weighted.high_pup_l.(cent_features{ni})-spon_states_weighted.low_pup_q.(cent_features{ni}),2));
%     figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%     plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%         [],  -L, L, 1,colormap(redblue))
%     title(['run vs pupil low ' cent_features{ni} ' Centrality (spont)']);
%     % mysave(gcf, fullfile(outputfiggolder,['spont_state3_minus_state1_',cent_features{ni},'_weighted_heatmap_allen']));
%     % trial
%     for state_i=1:length(statenames)
%         diffmask = scores_to_heatmap_allen(mean(trial_states_weighted.(statenames{state_i}).(cent_features{ni}).correct-trial_states_weighted.(statenames{state_i}).(cent_features{ni}).incorrect,2));
%         figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%         plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%             [],  -L, L, 1,colormap(redblue))
%         title(['Correct minus incorrect on ' statenames{state_i} ' ' cent_features{ni} ' Centrality (trial)']);
%
%         % mysave(gcf, fullfile(outputfiggolder,['trial_state' num2str(state_i) '_',cent_features{ni},'_weighted_heatmap_allen']));
%     end
%     % not weighted
%     %     graph_overlay_allen_paired('',fullfile('', 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%     %         spon_states_notweighted.low_pup_q.(cent_features{ni}),strcat('','/spon_run_pupillow/'),cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (spont)'],parcels_names,length(animals));
%     % spont
%     diffmask = scores_to_heatmap_allen(mean(spon_states_notweighted.high_pup_l.(cent_features{ni})-spon_states_notweighted.low_pup_q.(cent_features{ni}),2));
%     figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%     plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%         [],  -L, L, 1,colormap(redblue))
%     title(['run vs pupil low ' cent_features{ni} ' Centrality (spont)']);
%
%
%     % mysave(gcf, fullfile(outputfiggolder,['spont_state3_minus_state1_',cent_features{ni},'_notweighted_heatmap_allen']));
%     % trial not
%     for state_i=1:length(statenames)
%         diffmask = scores_to_heatmap_allen(mean(trial_states_notweighted.(statenames{state_i}).(cent_features{ni}).correct-trial_states_notweighted.(statenames{state_i}).(cent_features{ni}).incorrect,2));
%         figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%         plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%             [],  -L, L, 1,colormap(redblue))
%         title(['Correct minus incorrect on ' statenames{state_i} ' ' cent_features{ni} ' Centrality (trial)']);
%         % mysave(gcf, fullfile(outputfiggolder,['trial_state' num2str(state_i) '_',cent_features{ni},'_notweighted_heatmap_allen']));
%     end


% end

% end




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
% mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\','figure_1',name), 'all');
end


%%
function eval_diff_map(animals, statenames)
addpath(genpath('../centrality_measures/'));
outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality';
for ai=1:length(animals)
    animal=animals{ai};
    
    for state_i = 1:length(statenames)
        load(fullfile(outputfolder,'network_analysis_corr',statenames{state_i}),'W_corr');
        W_corr = threshold_cor_matrix(W_corr);
        W_corr = W_corr + eye(size(W_corr));
        [M(:, state_i, ai),Q(state_i, ai)]=community_louvain(abs(W_corr));
        P(:, state_i, ai)=participation_coef(abs(W_corr),M(:, state_i, ai),0);
        
        configParams.maxInd=20;
        [diffusion_map, Lambda, Psi, Ms, Phi, K_rw] = calcDiffusionMap(exp(W_corr), configParams);
        eigenvals(:, state_i, ai) = Lambda(2:end);
        firsteigenvec(:, state_i, ai) = Psi(:,2);
    end
end
M = mean(eigenvals, 3);
S = std(eigenvals, [],3)/sqrt(size(eigenvals,3)-1);
barwitherr(S,M)
end
% function plot_centrality_res_gal(animals, outputfiggolder, statenames)
% [~, ~, finalindex] = get_allen_meta_parcels;
% cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector', 'participation', 'community'};
% for state_i = 1:length(statenames)
%     for cent_i = 1:length(cent_features)
%         spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = zeros(256,256,length(animals));
%     end
% end
% for i=1:length(animals)
%     animal=animals{i};
%     outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality';
%     for state_i = 1:length(statenames)
%         load(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} 'LSSC']));
%         cent_features = fieldnames(cent_corr_notweighted);
%         for cent_i = 1:length(cent_features)
%             P = scores_to_heatmap_gal(cent_corr_notweighted.(cent_features{cent_i}), animal);
%             spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i})(:, :, i) = P;
%
%         end
%
%     end
% end
% mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
% for ni = 1:length(cent_features)-1
%     %difference maps, not weighted, for each centrality measure
%     braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
%     parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
%
%     diffmask = nanmean(spon_states_notweighted.high_pup_l.(cent_features{ni}),3)-...
%         nanmean(spon_states_notweighted.low_pup_q.(cent_features{ni}),3);
%     plot_heatmap(diffmask, -0.04, 0.04, ['Difference in Centrality ' cent_features{ni} 'high pup loc minus low pup q'], parcelsallen.parcells_new.indicators(:,:,finalindex));
%     % mysave(gcf, fullfile(outputfiggolder,'not_weighted','spon_run_pupillow',strcat(cent_features{ni},'_heatmap_gal')), 'all');
%     diffmask = nanmean(spon_states_notweighted.high_pup_q.(cent_features{ni}),3)-...
%         nanmean(spon_states_notweighted.low_pup_q.(cent_features{ni}),3);
%     plot_heatmap(diffmask, -0.04, 0.04, ['Difference in Centrality ' cent_features{ni} 'high pup q minus low pup q'], parcelsallen.parcells_new.indicators(:,:,finalindex));
%     % mysave(gcf, fullfile(outputfiggolder,'not_weighted','spon_pupilhigh_pupillow',strcat(cent_features{ni},'_heatmap_gal')), 'all');
% end
%
% end

% function plot_centrality_res(animals, outputfiggolder, statenames)
% [parcels_names] = get_allen_meta_parcels;
% cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector', 'participation', 'community'};
% for state_i = 1:length(statenames)
%     for cent_i = 1:length(cent_features)
%         spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = [];
%         %         spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = [];
%     end
% end
% for i=1:length(animals)
%     animal=animals{i};
%     outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality';
%     for state_i = 1:length(statenames)
%         load(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} 'Allen']));
%         cent_features = fieldnames(cent_corr_notweighted);
%         for cent_i = 1:length(cent_features)
%             spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
%                 spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}), ...
%                 cent_corr_notweighted.(cent_features{cent_i}));
%
%             %             spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
%             %                 spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}), ...
%             %                 cent_corr_weighted.(cent_features{cent_i}));
%         end
%
%     end
% end
% % mkNewDir(fullfile(outputfiggolder, 'weighted'))
% mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
% legstr = {'Low Q', 'High Q', 'Loc'};
% for ni = 1:length(cent_features)-1
%     graph_overlay_allen_2conditions([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.low_pup_q.(cent_features{ni}),...
%         spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%         'spon_threestates',cent_features{ni},['2 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr([1 3]));
%
%
%
%     %     graph_overlay_allen_3conditions([],fullfile(outputfiggolder, 'weighted'), spon_states_weighted.low_pup_q.(cent_features{ni}),...
%     %         spon_states_weighted.high_pup_q.(cent_features{ni}), spon_states_weighted.high_pup_l.(cent_features{ni}),...
%     %         'spon_threestates',cent_features{ni},['3 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr);
%     %
%     graph_overlay_allen_3conditions([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.low_pup_q.(cent_features{ni}),...
%         spon_states_notweighted.high_pup_q.(cent_features{ni}), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%         'spon_threestates',cent_features{ni},['3 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr);
%
%     %difference maps, not weighted, for each centrality measure
%     braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
%     parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
%
%     %high vs low pup
%     graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_q.(cent_features{ni}),...
%         spon_states_notweighted.low_pup_q.(cent_features{ni}),'spon_pupilhigh_pupillow',cent_features{ni},['pupil high vs low ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
%
%     graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
%         spon_states_notweighted.high_pup_q.(cent_features{ni}),spon_states_notweighted.low_pup_q.(cent_features{ni}),...
%         'spon_pupilhigh_pupillow',cent_features{ni},['pupil high vs low ' cent_features{ni} ' Centrality (spon)']);
%
%     %locomotion vs low pup
%     graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%         spon_states_notweighted.low_pup_q.(cent_features{ni}),'spon_run_pupillow',cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
%
%     graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
%         spon_states_notweighted.high_pup_l.(cent_features{ni}),spon_states_notweighted.low_pup_q.(cent_features{ni}),...
%         'spon_run_pupillow',cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (spon)']);
%
%     %locomotion vs high pup
%     graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%         spon_states_notweighted.high_pup_q.(cent_features{ni}),'spon_run_pupilhigh',cent_features{ni},['run vs pupil high ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
%
%     graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
%         spon_states_notweighted.high_pup_l.(cent_features{ni}),spon_states_notweighted.high_pup_q.(cent_features{ni}),...
%         'spon_run_pupilhigh',cent_features{ni},['run vs pupil high ' cent_features{ni} ' Centrality (spon)']);
%
% end
%
% end

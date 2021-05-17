function main_network_centrality_evaluateion_spont_crispr
addpath(genpath('../utils'));
addpath(genpath('D:/utils/affinity/'))
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
animals_db = get_animals_meta_data_by_csv;
statenames_3states = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
statenames_2states = {'qui', 'loc'};
statenames_5states = {'sit', 'loc','low_pupil','high_pupil','low_face','high_face'};

similarity_name = {'pearson_corr',  };%'corr',,  'L2' 'fullcorr' 'cov''partial_corr'
signames = { 'Allen'}; % 'Grid4'};
REPS = 1000;
for sim_i = 1:length(similarity_name)
    
    for ai = 1:length(animals_db.folder_list)
        if animals_db.toinclude_list(ai)==find(strcmp(animals_db.toinclude_lut, 'Good'))
%                      eval_weights_and_cent(signames, similarity_name{sim_i}, animals_db.folder_list{ai}, statenames_5states);
            eval_corrmat_and_shuffles(similarity_name{sim_i}, animals_db.folder_list{ai}, statenames_5states, REPS);
            %         eval_weights_and_cent_perday(similarity_name{sim_i}, animals{ai}, statenames);
            %         eval_weights_and_cent(signames, similarity_name{sim_i}, animals_db.folder_list{ai}, statenames_2states);
            %         eval_weights_and_cent(signames, similarity_name{sim_i}, animals_db.folder_list{ai}, statenames_3states);
            
        end
    end
end
end
%
% statenames = statenames_2states;
% for sig_i = 1:length(signames)
%     for sim_i = 1:length(similarity_name)
%         l=1;labels_type=[];labels_state=[];
%         outputfolder = ['X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr\network_centrality_' similarity_name{sim_i}];
%         for state_i = 1:length(statenames)
%
%             for ai = 1:length(animals_db.folder_list)
%                 if animals_db.toinclude_list(ai)==find(strcmp(animals_db.toinclude_lut, 'Good'))
%                     filename=fullfile(outputfolder,animals_db.folder_list{ai} ,[statenames{state_i} ,signames{sig_i} '_Inf.mat']);
%                     if ~isfile(filename)
%                         continue;
%                     end
%                    r= load(filename,'Wmat');
%                    for u=1:size(r.Wmat,3)
%                     Wvec(:,:,l) = r.Wmat(:,:,u)+eye(23);
%                     labels_type(l) = animals_db.type_list(ai);
%                     labels_state(l) = state_i;
%                     animal_label(l) = ai;
%                     cohort_label(l) = animals_db.cohort_list(ai);
%                     sex_label(l) = animals_db.sex_list(ai);
%                     l=l+1;
%                    end
%                 end
%             end
%         end
%
%         figure;i=1;
%         for sex_i = 2:3
%             for lab_i = 1:3
%                 selinds = setdiff(find(labels_type==lab_i&sex_label==sex_i&labels_state==2),1398)
%                 subplot(2,3,i);i=i+1;
%                 imagesc(mean(Wvec(:,:,selinds),3));
%                         title(['sex ' animals_db.sex_lut{sex_i} ' ' animals_db.type_lut{lab_i}]);
%
%             end
%         end
%         figure;i=1;
%         for sex_i = 2:3
%             selinds = find(sex_label==sex_i);
%
%             mRiemannianMean = RiemannianMean(Wvec(:,:,selinds));
%             mX = proj_R1(mRiemannianMean^(-1/2),Wvec(:,:,selinds));
%
%             par = SetGenericDimsQuestParams(2,false);
%             par.init_aff{2}.knn=max(round(size(mX,2)*0.01),2000);
%             [ initAll, sig ] = CalcInitAff2D( mX, par.init_aff{2} );
%             %     [ initAll, sig ] = CalcInitAff2D( mX(:,cohort_label==3&sex_label==2), par.init_aff{2} );
%             dParams.maxInd = min(size(initAll,1),100);
%             dParams.normalization='lb';
%             [diffmap_long_p, Lambda] = calcDiffusionMap(initAll,dParams);
%             figure;subplot(2,2,1);plotEmbeddingWithColors(diffmap_long_p(1:3,:)',labels_type(selinds));title('Type');
%             subplot(2,2,2);plotEmbeddingWithColors(diffmap_long_p(1:3,:)',labels_state(selinds));title('state');
%             subplot(2,2,3);plotEmbeddingWithColors(diffmap_long_p(1:3,:)',animal_label(selinds));title('animal');
%             suptitle(['sex ' num2str(sex_i) ' co ' num2str(lab_i)]);
%
%         end
%         %
% %
% %
% %     figure;subplot(2,2,1);plotEmbeddingWithColors(diffmap_long_p',labels_type);
% %     i=labels_state==1;
% %     subplot(2,2,2);plotEmbeddingWithColors(diffmap_long_p(:,i)',labels_type(i));
% %     i=labels_state==2;
% %     subplot(2,2,3);plotEmbeddingWithColors(diffmap_long_p(:,i)',labels_type(i));
% %
% %     figure;subplot(2,2,1);plotEmbeddingWithColors(diffmap_long_p',labels_state);
% %     i=labels_type==1;
% %     subplot(2,2,2);plotEmbeddingWithColors(diffmap_long_p(:,i)',labels_state(i));
% %     i=labels_type==2;
% %     subplot(2,2,3);plotEmbeddingWithColors(diffmap_long_p(:,i)',labels_state(i));
% %      i=labels_type==3;
% %     subplot(2,2,4);plotEmbeddingWithColors(diffmap_long_p(:,i)',labels_state(i));
% %
%     end
% end
% end

function eval_weights_and_cent_perday(simname, animal, statenames)
outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality_' simname];


mkNewDir(outputfolder);

[~, allen_parcels] = getParcellsByLansAllansAtlas;
[parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;

regionLabel.Allen = allen_parcels.regionNum;
regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
% [roiLabelsbyAllen_gal, regionLabel.Gal, maskByAllen_gal, maskByAllen.Gal] = get_gal_parcels_lables(animal);
[parcels_names.LSSC, ~, ~, regionLabel.LSSC] = get_gal_meta_parcels_by_allen(parcels_names.Allen, maskByAllen.Allen, ...
    regionLabel.Allen, animal);

load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\',animal,'_spont_data_3states_dfff.mat'],...
    'low_pup_q','high_pup_q','high_pup_l'); %#ok<NASGU>
disp(animal);
signames = {'Grid4','Allen'};
[~,days2process] = animaltodays(animal);
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    if 0&&exist(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,'Allen.mat']),'file') &&...
            exist(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,'LSSC.mat']),'file')
        continue;
    end
    
    data = eval(statenames{state_i});
    for sig_i = 1:length(signames)
        for dd=1:length(days2process)
            disp([animal ' ' signames{sig_i} ' ' num2str(days2process(dd)) ' ' statenames{state_i}]);
            currdata = data.(signames{sig_i});
            currdata=currdata(:, data.days==days2process(dd));
            if sum(data.days==days2process(dd))==0
                continue;
            end
            if any(currdata(:))
                tt = ~isnan(sum(currdata));
                if sum(tt)==0
                    disp('nans in dataset')
                    
                    continue;
                else
                    currdata = currdata(:,tt);
                end
            end
            if exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_day' num2str(days2process(dd)) '.mat']),'file')
                load(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_day' num2str(days2process(dd)) '.mat']),'W_corr')
            else
                W_corr = measure_weights(currdata, simname);
            end
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}));
            save(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signames{sig_i} '_day' num2str(days2process(dd)) '.mat']),'W_corr',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
            
        end
    end
end
end
function eval_corrmat_and_shuffles(simname, animal, statenames, REPS)

outputfolder = ['X:\Hadas\Meso-imaging\CRISPR\analysis_results\network_centrality_' simname];

mkNewDir(outputfolder);

[~, allen_parcels] = getParcellsByLansAllansAtlas;
[parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;

regionLabel.Allen = allen_parcels.regionNum;
regionLabel.Allen=regionLabel.Allen(finalindex.Allen);


if ~exist(fullfile('X:\Hadas\Meso-imaging\CRISPR\analysis_results\',animal,'con_states.mat'), 'file')
    return;
end
disp(animal);
mkNewDir(fullfile(outputfolder,animal));

res1 = load(fullfile('X:\Hadas\Meso-imaging\CRISPR\analysis_results\',animal,'arousal_traces_states.mat'));

xa = load(fullfile('X:\Hadas\Meso-imaging\CRISPR\traces_data', animal, 'Ca_traces_spt_patch11_Allen_dfff'));
parcels_traces=xa.parcels_time_trace(finalindex.Allen, :);

load(fullfile('X:\Hadas\Meso-imaging\CRISPR\traces_data', animal, 'smrx_signals_v4.mat'),'timing');
t_imaging = timing.bluestart;
L = min([length(t_imaging) size(parcels_traces,2)]);
parcels_traces=parcels_traces(:,1:L);
t_imaging=t_imaging(1:L);
notnans = all(~isnan(parcels_traces));
parcels_traces=parcels_traces(:,notnans);

t_imaging = t_imaging(notnans);
corrmat_highpupilloc=[];corrmat_pupil=[];

if isfield(res1.segments_arousals,'loc') && isfield(res1.segments_arousals,'high_pupil')&&...
        ~isempty(res1.segments_arousals.loc)&&~isempty(res1.segments_arousals.sit)
    corrmat_highpupilloc = extract_shuffled_segs(simname, res1.segments_arousals, t_imaging, parcels_traces, 'loc', 'high_pupil', REPS);
end
if isfield(res1.segments_arousals,'low_pupil') && isfield(res1.segments_arousals,'high_pupil')&&...
        ~isempty(res1.segments_arousals.low_pupil)&&~isempty(res1.segments_arousals.high_pupil)
    corrmat_pupil = extract_shuffled_segs(simname, res1.segments_arousals, t_imaging, parcels_traces, 'high_pupil', 'low_pupil', REPS);
end


save(fullfile(outputfolder,animal ,['shuffled_corr_Allen_.mat']),...
    'corrmat_highpupilloc','corrmat_pupil');
end


function eval_weights_and_cent(signames, simname, animal, statenames)
outputfolder = ['X:\Hadas\Meso-imaging\CRISPR\analysis_results\network_centrality_' simname];

datafilename = ['spont_data_' num2str(length(statenames)) 'states_dfff.mat'];

mkNewDir(outputfolder);
for s=1:length(signames)
    switch signames{s}
        case 'Allen'
            [~, allen_parcels] = getParcellsByLansAllansAtlas;
            [parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
            
            regionLabel.Allen = allen_parcels.regionNum;
            regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
        case 'LSSC'
            % [roiLabelsbyAllen_gal, regionLabel.Gal, maskByAllen_gal, maskByAllen.Gal] = get_gal_parcels_lables(animal);
            [parcels_names.LSSC, ~, ~, regionLabel.LSSC] = get_gal_meta_parcels_by_allen(parcels_names.Allen, maskByAllen.Allen, ...
                regionLabel.Allen, animal);
        case 'Grid4'
            load('X:\Hadas\Meso-imaging\lan\xspsych\spt\xs_31_grid4_dfff.mat','par_inds');
            [parcels_names.Grid4, regionLabel.Grid4, finalindex.Grid4, regionLabel.nameslegend, maskByAllen.Grid4, labelsbyallen.Grid4] = getAllenClusteringLabelsGrid(par_inds, 4);
            ii = discard_inds;
            
            parcels_names.Grid4=parcels_names.Grid4(ii);
            regionLabel.Grid4=regionLabel.Grid4(ii);
    end
end
if ~exist(fullfile('X:\Hadas\Meso-imaging\CRISPR\analysis_results\',animal,'con_states.mat'), 'file')
    return;
end
load(fullfile('X:\Hadas\Meso-imaging\CRISPR\analysis_results\',animal,'con_states.mat')); %#ok<NASGU>
disp(animal);
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    mkNewDir(fullfile(outputfolder,animal));
    data = eval(statenames{state_i});
    if isempty(data.t)
        continue;
    end
    for sig_i = 1:length(signames)
        
        data.(signames{sig_i}) = data.(signames{sig_i})(:, all(~isnan(data.(signames{sig_i}))));
        if any(isnan(data.(signames{sig_i})(:)))
            tt = ~isnan(sum(data.(signames{sig_i})));
            if sum(tt)==0
                disp('nans in dataset')
                
                continue;
            else
                data.(signames{sig_i}) = data.(signames{sig_i})(:,tt);
            end
        end
        
        
        if strcmp(signames{sig_i}, 'Grid4')
            if length(ii) > size(data.(signames{sig_i}),1)
                error('check this');
            else
                data.(signames{sig_i}) = data.(signames{sig_i})(ii,:);
            end
            thT=[Inf ];
        else
            thT=[Inf  ];
        end
        %         if  exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '.mat']),'file')
        %             load(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '.mat']),'W_corr')
        %         else
        [W_corr,Wmat] = measure_weights_bysegs(1:length(data.t),data.(signames{sig_i}), simname);
        if isempty(W_corr)
            continue;
        end
        %         end
        for th=thT
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}), @process_sim, th);
            
            %         [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}));
            save(fullfile(outputfolder,animal ,[statenames{state_i} ,signames{sig_i} '_' num2str(th) '.mat']),'W_corr',...
                'cent_corr_weighted','Wmat',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
        end
    end
    
end
end






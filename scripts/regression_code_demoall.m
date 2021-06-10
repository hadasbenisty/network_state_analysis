function regression_code_demoall

% cd('C:\Users\SuperComp-HigleyLab\Documents\Lav\network_state_analysis\scripts')
addpath(genpath('../pre_processing'));
addpath(genpath('../meta_data_processing'));
addpath(genpath('../network_state_analysis/utils'));
addpath(genpath('../network_state_analysis/functions'));
addpath(genpath('../network_state_analysis/scripts'));

datapath = 'X:\Hadas\Meso-imaging\CRISPR\traces_data';
imagingpath='X:\Hadas\Meso-imaging\CRISPR\traces_data';
spike2path = 'X:\Hadas\Meso-imaging\CRISPR\traces_data';
animals_db = get_animals_meta_data_by_csv;
%procdatapath = 'X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr\regression';
procdatapath='X:\Hadas\Meso-imaging\CRISPR\analysis_results';
%figspath = 'X:\Hadas\Mesoimaging\crispr\meso_results\figs_crispr\regression';
figspath='X:\Hadas\Meso-imaging\CRISPR\Figures\regression';
do_analysis(animals_db, datapath, imagingpath, spike2path, procdatapath);
plot_regression_results(animals_db, procdatapath, figspath);

end
function do_analysis(animals_db, datapath, imagingpath, spike2path, procdatapath)
% regression analysis
addpath(genpath('../network_state_analysis/'));
clc;
for k=1:length(animals_db.animal_list)
    if animals_db.toinclude_list(k)==find(strcmp(animals_db.toinclude_lut, 'Good'))&&...
            animals_db.isgoodpupil_list(k)==find(strcmp(animals_db.isgoodpupil_lut, 'GOOD'))
        mkNewDir(fullfile(procdatapath,animals_db.folder_list{k}));
        resfile = fullfile(procdatapath,animals_db.folder_list{k}, 'regression_res.mat');
        if isfile(resfile)&&0
            a = load(resfile);
            if ~isempty(a.outputMat_gcamp_regular)
                continue;
            end
        end
        %extract_sustained_state(datapath, procdatapath, spike2path, animals_db.folder_list{k});
        tic
        isair=animals_db.sessionsid_list(k)==1;
%         try
            outputMat_gcamp_regular = ...
                returnCVRSquared(animals_db.folder_list{k},k,fullfile(datapath,animals_db.folder_list{k}),fullfile(imagingpath,animals_db.folder_list{k}),fullfile(spike2path,animals_db.folder_list{k}),isair,25);
%         catch
%             continue;
%         end
        toc
        if ~isempty(outputMat_gcamp_regular)
            save(resfile,'outputMat_gcamp_regular');
        end
        
    else
        %         disp(strcat('data isnt good for',animals_db.folder_list{k}));
    end
    
end
end
function plot_regression_results(animals_db, procdatapath, figspath)
parcels_names = get_allen_meta_parcels;
all_results_type = nan(4, length(parcels_names),  length(animals_db.folder_list));
for k=1:length(animals_db.folder_list)
    resfile = fullfile(procdatapath,animals_db.folder_list{k}, 'regression_res.mat');
    if ~isfile(resfile)
        continue;
    end
    load(resfile,'outputMat_gcamp_regular');
    if isempty(outputMat_gcamp_regular)
        continue;
    end
    all_results_type(:,:,k) = outputMat_gcamp_regular;
   
    
end



M=0.35;
for ti = 1:length(animals_db.type_lut)
    curtype = animals_db.type_lut{ti};
    %order: allVars,pupil,face,wheel;
    %all
    rsquareRes = all_results_type(:, :, animals_db.type_list==ti);
    
    sessionsnum=sum(~isnan(squeeze(rsquareRes(1,1,:))));
    %% plot all sessions and animals
    R = ceil(sqrt(sessionsnum+1));
    Pvec = scores_to_heatmap_allen(squeeze(rsquareRes(1,:,:)), false);
    figure;inds=find(~isnan(squeeze(rsquareRes(1,1,:))));
    animalsinds = find(animals_db.type_list==ti);
    for j=1:sessionsnum
       subplot(R,R,j);imshow(Pvec(:,120:240,inds(j)),[0 M]);colorbar; 
       ttl = animals_db.folder_list{animalsinds(j)};
       ttl = ttl(strfind(ttl,'/')+1:end);
       ttl(ttl=='_') = ' ';
       title( ttl);
    end
    subplot(R,R,j+1);imshow(nanmean(Pvec(:, 120:240,:),3),[0 M]);title('Average');colorbar;colormap jet;
     suptitle([ curtype '  per session, using all predictors']);
    figure;subplot(2,2,1);
    sessionsnum=sum(~isnan(squeeze(rsquareRes(1,1,:))));
    imagesc(nanmean(Pvec(:, 120:240,:),3),[0 M]);title(['All N=' num2str(sessionsnum)]);colorbar;
    
    %Pupil
    subplot(2,2,2);
    Pvec = scores_to_heatmap_allen(squeeze(rsquareRes(2,:,:)), false);
    sessionsnum=sum(~isnan(squeeze(rsquareRes(2,1,:))));
    imagesc(nanmean(Pvec(:, 120:240,:),3),[0 M]);title(['Pupil N=' num2str(sessionsnum)]);colorbar;
    
    %Face
    subplot(2,2,3);
    Pvec = scores_to_heatmap_allen(squeeze(rsquareRes(3,:,:)), false);
    sessionsnum=sum(~isnan(squeeze(rsquareRes(3,1,:))));
    imagesc(nanmean(Pvec(:, 120:240,:),3),[0 M]);title(['Face N=' num2str(sessionsnum)]);colorbar;
    
    %Wheel
    subplot(2,2,4);
    Pvec = scores_to_heatmap_allen(squeeze(rsquareRes(4,:,:)), false);
     sessionsnum=sum(~isnan(squeeze(rsquareRes(4,1,:))));
    imagesc(nanmean(Pvec(:, 120:240,:),3),[0 M]);title(['Wheel N=' num2str(sessionsnum)]);colorbar;
   
    
    colormap jet;
    suptitle([ curtype ' ' num2str(sessionsnum),' sessions']);
    
    mysave(gcf, fullfile(figspath, ['R2_' curtype]));
end
end
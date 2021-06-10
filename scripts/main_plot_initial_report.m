function main_plot_initial_report
%
addpath('..\functions');
addpath(genpath('../utils'));
addpath(genpath('../meta_data_processing/'));
% procdatapath = 'X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr';
[~, ~, finalindex.Allen] = get_allen_meta_parcels;
outputfiggolder = 'X:\Hadas\Meso-imaging\CRISPR\Figures\initial_report';
dffpath = 'X:\Hadas\Meso-imaging\Antara\final_preprocess\alldata';
dffpath = 'X:\Hadas\Meso-imaging\CRISPR\traces_data\';

mkNewDir(outputfiggolder);
animals_db = get_animals_meta_data_by_csv;
for fi = 1:length(animals_db.ex_lut)
    mkNewDir(fullfile(outputfiggolder , ['group' num2str(animals_db.ex_lut(fi))]));
end
parcelsallen=load('parcells_updated121519.mat');

spontindic = animals_db.sessionsid_list==find(strcmp(animals_db.sessionsids_lut, 'Spon'));
for i = 1:length(animals_db.animal_lut)
    curranimal = find((animals_db.animal_list ==  (i))&spontindic);
    animals_inds_list{i} = curranimal;
end
for i = 1:length(animals_inds_list)
    
    for ii = 1:length(animals_inds_list{i})
        k = animals_inds_list{i}(ii);
        if  animals_db.sessionsid_list(k)==find(strcmp(animals_db.sessionsids_lut, 'Air'))
            continue;
        end
        outfigscurrpath = fullfile(outputfiggolder , ['group' num2str(animals_db.ex_lut(animals_db.ex_list(k)))]);
        
        str = animals_db.animal_lut{animals_db.animal_list(k)};
        str(str == '_') = '-';
        if isfile(fullfile(outfigscurrpath, [str '_means.fig']))
            continue;
        end
% spike2pth = fullfile('X:\Hadas\Meso-imaging\Antara\data\Antara\AnalyzedData\', animals_db.folder_list{k});
        spike2pth = fullfile('X:\Hadas\Meso-imaging\CRISPR\traces_data\');
        if ~isfile(fullfile(spike2pth, animals_db.folder_list{k},'smrx_signals_v4.mat'))
            disp(['No spike2 '  animals_db.folder_list{k}]);
            continue;
        end
        datafile_allen = fullfile(dffpath, animals_db.folder_list{k},  'Ca_traces_spt_patch11_Allen_dfff.mat');
        if ~isfile(datafile_allen)
            disp(['No allen '  animals_db.folder_list{k}]);
            continue;
        end
        fulldfof = fullfile(dffpath, animals_db.folder_list{k}, 'Ca_traces_spt_patch11_dfff.mat');
        if ~isfile(fulldfof)
            disp(['No spt '  animals_db.folder_list{k}]);
            continue;
        end
        rawblue = fullfile(dffpath, animals_db.folder_list{k}, 'detrended_blue_dfff_al.mat');
        if ~isfile(rawblue)
            disp(['No raw '  animals_db.folder_list{k}]);
            
            continue;
        end
        %% plot traces

        load(fullfile(spike2pth,  animals_db.folder_list{k},'smrx_signals_v4.mat'), 'timing', 'channels_data');
        load(datafile_allen, 'parcels_time_trace');
        t_imaging = timing.bluestart;
        L=min(length(t_imaging), size(parcels_time_trace,2));
        figure;
        h(1) = subplot(3,1,1);plot((1:length(channels_data.wheelspeed))/5000, channels_data.wheelspeed)
        title('Wheel');
        h(2) = subplot(3,1,2);plot(t_imaging(1:L), parcels_time_trace(2,1:L));title('V1');
        hold all;plot(t_imaging(1:L), parcels_time_trace(1,1:L));legend('left','right');
        ylabel('\Delta f/f');
        h(3) = subplot(3,1,3);plot(t_imaging(1:L), parcels_time_trace(52,1:L));title('M1');ylabel('\Delta f/f');
        hold all;plot(t_imaging(1:L), parcels_time_trace(51,1:L));legend('left','right');
        linkaxes(h,'x');
        xlabel('Time [samples]');
        str = animals_db.animal_lut{animals_db.animal_list(k)};
        str(str == '_') = '-';suptitle(str);
        try
            mysave(gcf, fullfile(outfigscurrpath, [str '__traces']));
        catch
            mysave(gcf, fullfile(outfigscurrpath, [str '___traces.jpg']));
        end
        corr_right = corr(parcels_time_trace(finalindex.Allen,1e2:end-100)');
        corr_left = corr(parcels_time_trace(finalindex.Allen-1,1e2:end-100)');
        
        load('brain_mask');
        load(fulldfof, 'regsig');
        load(rawblue, 'dff_blue','mval','sval');
        P=zeros(256);P(brain_mask==1) = dff_blue(:,1000);
        figure;subplot(3,3,3);imagesc(P);title('zscored blue frame');colorbar;
        P(brain_mask==1) = mval;
        subplot(3,3,1);imagesc(P);title('mean raw blue blue');colorbar;
        P(brain_mask==1) = sval;
        subplot(3,3,2);imagesc(P);title('std raw blue pixels');colorbar;
        subplot(3,3,4);imagesc(reshape(regsig(:,1e3),[256 256]));title('regressed frame');colorbar;
        
        subplot(3,1,3);bar(nanstd(parcels_time_trace,[],2));title('std parcels');
        set(gca,'XTick',1:56);set(gca,'XTickLabel',parcelsallen.parcells_new.names);
        subplot(3,3,5);imagesc(corr_right,[0 1]);colorbar;title('Parcels corr right');
        subplot(3,3,6);imagesc(corr_left,[0 1]);colorbar;title('Parcels corr left');
        str = animals_db.animal_lut{animals_db.animal_list(k)};
        str(str == '_') = '-';suptitle(str);
        set(gcf,'Position',[1          41        1920         963]);
        mysave(gcf, fullfile(outfigscurrpath, [str '_means']));
        close all;
        break;
    end
end
end
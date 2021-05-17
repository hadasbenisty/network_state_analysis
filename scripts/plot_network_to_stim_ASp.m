function plot_network_to_stim_ASp
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../utils/'))
addpath(genpath('../../utils/Questionnaire/'))
addpath(genpath('../network_analysis/'));
figsfolder = 'X:\Hadas\Meso-imaging\GRABS_results\Figures\network_analysis';
proc_output_folder = 'X:\Hadas\Meso-imaging\GRABS_results\processing_dir\';
inputFolder='X:\GRABS_Data\Analyzed_SVDMethodPatch14\DualMice';
clc;
MiceAnalyze=5:10;
winLenvec=[1 3 5 10 15];
names = {'Allen_green','Allen_blue'};
Acc_phi = nan(99,6, length(winLenvec), length(MiceAnalyze), 8, 2);
Acc_im = nan(6, length(winLenvec), length(MiceAnalyze), 8, 2);
Acc_im_parcels = nan(23, 6, length(winLenvec), length(MiceAnalyze), 8, 2); 
for ni=length(names):-1:1
    for animal_ind=1:6
        animal = MiceAnalyze(animal_ind);
        mainDir1 =fullfile(inputFolder, sprintf("grabAM%02d",(animal)), 'imaging with 575 excitation');
        folders = dir(mainDir1);
        folders=folders(3:end);
        for folder = 1:length(folders)
            for wi=1:length(winLenvec)
                
                
                
                currfolder = fullfile(mainDir1, folders(folder).name);
                if  contains(currfolder, 'PostDrug')||~contains(currfolder,'vis')
                    continue;
                end
                
                
                
                
                
                accresfile=fullfile(proc_output_folder, 'network_analysis',['acc_' names{ni} '_' folders(folder).name '_sims' num2str(winLenvec(wi)) '.mat']);
                if ~isfile(accresfile)
                    disp(accresfile);
                    continue;
                end
                
                
                
                load(accresfile, 'accuracy_phi','accuracy_im','accuracy_im_parcels');
                Acc_phi(:, :, wi, animal_ind, folder, ni) = accuracy_phi';
                Acc_im(:, wi, animal_ind, folder, ni) = accuracy_im;
                Acc_im_parcels(:,  :, wi, animal_ind, folder, ni) = accuracy_im_parcels';
            end
            
        end
        
    end
end
%%
ttl={'Activity','\phi_c'};
clrs = ['g','b'];
N=6;
xtmp = squeeze(nanmean(nanmean(Acc_im,2),4));

figure;subplot(2,1,1);
for ni=1%:2
    
Acc_im_M = nanmean(xtmp(:,:,ni),2);
Acc_im_S = nanstd(xtmp(:,:,ni),[],2)/sqrt(N);
shadedErrorBar([2 5 10 20 50 100],Acc_im_M, Acc_im_S,'lineprops',clrs(ni));
hold all
end

for ni=1%:2
    xtmp = nanmean(Acc_phi(:,:,:,:,:,ni),5);

Acc_phi_M = permute(nanmean(xtmp,4),[2 3 1]);
Acc_phi_S = permute(nanstd(xtmp,[],4)/sqrt(N),[2 3 1]);


shadedErrorBar([2 5 10 20 50 100],Acc_phi_M(:,2,20), Acc_phi_S(:,2,20),'lineprops',clrs(2));
hold all;
end

set(gca, 'XScale','log');
xlim([2 100]);ylabel('Accuracy')
ylim([0 1]);xlabel('Contrast');set(gca,'YTick',[0 0.5 0.6 1])
line(get(gca,'XLim'), [1 1]*0.5, 'LineStyle','--','Color','k')
c=get(gca,'Children');legend(c(end:-4:1),{'activity','\phi_c'},'Location','BestOutside');
L=[.45 .55];
load('rb');
xtmp = squeeze(nanmean(nanmean(Acc_im_parcels,5),3));
for ni=1%:2

    
    
Acc_im_M = nanmean(xtmp(:,:,:,ni),3);
subplot(2,2,3);
Pvec = scores_to_heatmap_allen(Acc_im_M, 0);
 plot_vals_heatmap(Pvec(:,:,end-1),'',[],  L(1), L(2), 1,rb);
            
colorbar;
title(['Activity Per Parcels ' ttl{ni}]);
 
end
%mysave(gcf,fullfile(figsfolder, 'stim_prediction'));
%%


for ni=1:2
    xtmp = nanmean(Acc_phi(:,:,:,:,:,ni),5);

Acc_phi_M = permute(nanmean(xtmp,4),[2 3 1]);
Acc_phi_S = permute(nanstd(xtmp,[],4)/sqrt(N),[2 3 1]);
subplot(3,1,3);

shadedErrorBar([2 5 10 20 50 100],Acc_phi_M(:,2,20), Acc_phi_S(:,2,20),'lineprops',clrs(ni));
hold all;
end
title('Correlation - 3 seconds window');
set(gca, 'XScale','log');
xlim([2 100]);ylim([0 1]);
line(get(gca,'XLim'), [1 1]*0.5, 'LineStyle','--','Color','k')
c=get(gca,'Children');legend(c(end:-4:1),ttl,'Location','BestOutside');
mysave(gcf,fullfile(figsfolder, 'stim_prediction1'));

%%
ttl={'RCaMP','AC'};
clrs = ['g','b'];
N=6;
xtmp = squeeze(nanmean(nanmean(Acc_im,2),4));

figure;subplot(3,1,1);
for ni=1%:2
    
Acc_im_M = nanmean(xtmp(:,:,ni),2);
Acc_im_S = nanstd(xtmp(:,:,ni),[],2)/sqrt(N);
shadedErrorBar([2 5 10 20 50 100],Acc_im_M, Acc_im_S,'lineprops',clrs(ni));
hold all
end
set(gca, 'XScale','log');
xlim([2 100]);ylabel('Accuracy')
ylim([0 1]);xlabel('Contrast');set(gca,'YTick',[0 0.5 0.6 1])
line(get(gca,'XLim'), [1 1]*0.5, 'LineStyle','--','Color','k')
title('All Parcels');c=get(gca,'Children');legend(c(end:-4:1),ttl,'Location','BestOutside');
L=[.45 .55];
load('rb');
xtmp = squeeze(nanmean(nanmean(Acc_im_parcels,5),3));
for ni=1:2

    
    
Acc_im_M = nanmean(xtmp(:,:,:,ni),3);
subplot(3,2,ni+2);
Pvec = scores_to_heatmap_allen(Acc_im_M, 0);
 plot_vals_heatmap(Pvec(:,:,end-1),'',[],  L(1), L(2), 1,rb);
            
colorbar;
title(['Per Parcels ' ttl{ni}]);
 
end
%mysave(gcf,fullfile(figsfolder, 'stim_prediction'));
%%


for ni=1:2
    xtmp = nanmean(Acc_phi(:,:,:,:,:,ni),5);

Acc_phi_M = permute(nanmean(xtmp,4),[2 3 1]);
Acc_phi_S = permute(nanstd(xtmp,[],4)/sqrt(N),[2 3 1]);
subplot(3,1,3);

shadedErrorBar([2 5 10 20 50 100],Acc_phi_M(:,2,20), Acc_phi_S(:,2,20),'lineprops',clrs(ni));
hold all;
end
title('Correlation - 3 seconds window');
set(gca, 'XScale','log');
xlim([2 100]);ylim([0 1]);
line(get(gca,'XLim'), [1 1]*0.5, 'LineStyle','--','Color','k')
c=get(gca,'Children');legend(c(end:-4:1),ttl,'Location','BestOutside');
mysave(gcf,fullfile(figsfolder, 'stim_prediction1'));

% subplot(4,1,4);
% 
% for ni=1:2
%     xtmp = nanmean(Acc_phi(:,:,:,:,:,ni),5);
% 
% Acc_phi_M = permute(nanmean(xtmp,4),[2 3 1]);
% Acc_phi_S = permute(nanstd(xtmp,[],4)/sqrt(N),[2 3 1]);
% 
% 
% shadedErrorBar([2 5 10 20 50 100],Acc_phi_M(:,3,20), Acc_phi_S(:,3,20),'lineprops',clrs(ni));
% hold all;
% end
% set(gca, 'XScale','log');
% xlim([2 100]);ylim([0 1]);
% xlabel('Contrast');
% 
% 
% line(get(gca,'XLim'), [1 1]*0.5, 'LineStyle','--','Color','k')
% legend(c(end:-4:1),ttl);
% title('Correlation - 5 seconds window');
% c=get(gca,'Children');legend(c(end:-4:1),ttl,'Location','BestOutside');
% 
% disp('finished');

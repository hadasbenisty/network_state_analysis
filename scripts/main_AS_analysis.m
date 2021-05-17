addpath(genpath('../gspbox/'));
addpath(genpath('../centrality_measures/'));
addpath(genpath('../meta_data_processing/AS/'));
addpath(genpath('../utils/'))
% clear;
addpath('..\meta_data_processing');
figuresFolder='X:\Hadas\Meso-imaging\GRABS_results\Figures\StateSpontCorrelationsLeftHemPermutationTest';
proc_output_folder = 'X:\Hadas\Meso-imaging\GRABS_results\processing_dir\';
inputFolder='X:\GRABS_Data\Analyzed_SVDMethodPatch14\DualMice';
MiceAnalyze=[{'grabAM05\imaging with 575 excitation\'},{'grabAM06\imaging with 575 excitation\'},{'grabAM07\imaging with 575 excitation\'},{'grabAM08\imaging with 575 excitation\'},{'grabAM09\imaging with 575 excitation\'},{'grabAM10\imaging with 575 excitation\'}];
[parcels_names, parcels_region_labels, finalindex, region_lut, allen_map_final_index, parcels_namesall] = get_allen_meta_parcels;
imaging_data.region_lut = region_lut;
Condition='NoDrug'; %'NoDrug','PreDrug','PostDrug' ol'}]; % whethere it's a drug or drug free session;
Np = length(parcels_names);

%get the size of each parcell so we can scale each cell in the color map accrdingly
pixelSizeParcell = hist(allen_map_final_index(:),unique(allen_map_final_index));
pixelSizeParcell=pixelSizeParcell(2:end);
totalBrainPixels=sum(pixelSizeParcell);
pixelpropSizeParcell=pixelSizeParcell/totalBrainPixels;
V1Idx=1; S1bIdx=14; M2Idx=23; %updated parcell idx for left V1, S1, and M2
load(fullfile(proc_output_folder,['IndivMouseOutput' Condition '.mat']),'allFaceLow','allFaceHigh',...
    'allPupHigh','allPupLow','allLoc','allSit',...
    'PupilHighImaging','PupilLowImaging','FaceHighImaging','FaceLowImaging',...
    'locImaging','sitImaging','CompData');

%% for each mouse, load spont and airpuff folders and perform correlations
Wloc = nan(Np, Np, size(PupilLowImaging, 1), size(locImaging, 2));
Wsit = nan(Np, Np, size(PupilLowImaging, 1), size(locImaging, 2));
networknums = [1 2 5 6 7];
imaging_data.fsample = 10;%Hz
winsizeSec = 8;
winhopSec = 4;
for animal=1:size(PupilLowImaging, 1)
    
    %for each subfolder
    for folder =1:size(PupilLowImaging, 2)
               %% time traces
        mainDir1 =fullfile(inputFolder,MiceAnalyze{animal});
        [dFoF_parcells, imaging_data.time, behaveData.wheel, behaveData.pupil, behaveData.facemap, ...
    behaveData.wheel_time, behaveData.pupil_time, LSSC_regionLabelsAllen, LSSC_roiLabelsByAllen, LSSC_parcelsMask] = get_AS_data(mainDir1, folder, Condition);

%% Allen
imaging_data.data = dFoF_parcells.Allen;
imaging_data.regionLabelsAllen = parcels_region_labels;
imaging_data.roiLabelsByAllen = 1:length(parcels_names);
imaging_data.parcelsMask = allen_map_final_index;
winsizeSec=30;%[3 5 10 15 20 25 30 ];
winhopSec=20;%1.5;
for wi=1:length(winsizeSec)
[results_allen_allendata{wi}, results_tree_allendata{wi}] = network_analysis_single_session(imaging_data, behaveData, winsizeSec(wi), winhopSec(wi),false);
end

%% LSSC
imaging_data.data = dFoF_parcells.LSCC;
imaging_data.regionLabelsAllen = LSSC_regionLabelsAllen;
imaging_data.roiLabelsByAllen = LSSC_roiLabelsByAllen;
imaging_data.parcelsMask = LSSC_parcelsMask;
winsizeSec=[3 5 10 15 20 25 30 ];
winhopSec=1.5;
for wi=1:length(winsizeSec)

[results_allen_LSSCdata{wi}, results_tree_LSSCdata{wi}] = network_analysis_single_session(imaging_data, behaveData, winsizeSec(wi), winhopSec(wi), false);
end

        pupil_down = interp1(pupil_time,pupil_Norm, imaging_time(t_win));
        wheel_down = interp1(wheel_time, wheel_speed, imaging_time(t_win));
        cc_p_tv = corr( pupil_down, TV');
        [~,ic]=sort(abs(cc_p_tv), 'descend');
        cc_p_corsum = corr( pupil_down, corr_sum');
        cc_p_tv_perm=[];cc_p_corsum_perm=[];
        for r=1:10
            cc_p_tv_perm(:,r) = corr( pupil_down(randperm(length(pupil_down))), TV');
            cc_p_corsum_perm(:,r) = corr( pupil_down(randperm(length(pupil_down))), corr_sum');
        end
        cc_p_tv_perm=mean(cc_p_tv_perm,2)';
        cc_p_corsum_perm=mean(cc_p_corsum_perm,2)';
        figure;bar(abs([cc_p_tv(ic);cc_p_corsum(ic);cc_p_tv_perm(ic);cc_p_corsum_perm(ic)]'))
        legend('tv','sum(corr)','perm tv','perm sum(corr)');ylabel('ABS Correlation with Pupil');
        set(gca,'XTick',1:length(brain_areas));set(gca,'XTickLabel',brain_areas)
        ylim([0 0.5]);
        cenmat=cenmat(~strcmp(centralitynames, 'pagerank')&~strcmp(centralitynames,'diffmap'),:,:);
        centralitynames=centralitynames(~strcmp(centralitynames, 'pagerank')&~strcmp(centralitynames,'diffmap'));
        for kk=1:size(cenmat,1)
            cc_p_cent(kk, :) = corr(pupil_down, squeeze(cenmat(kk,networknums,:))');
        end
        figure;bar(abs(cc_p_cent'))
        legend(centralitynames);ylabel('ABS Correlation with Pupil');
        set(gca,'XTick',1:length(networknums));set(gca,'XTickLabel',region_lut(networknums))
        ylim([0 0.5]);
        for ni=1:length(centralitynames)
            figure;
            wheel_down_s = filter(ones(winsize,1)/winsize, 1, wheel_down);
            wheel_down_s=[wheel_down_s(winsize/2+1:end); zeros(winsize/2,1)];
            hh(1) = subplot(3,1,1);
            plot(wheel_time,wheel_speed);
            hold all;
            plot(imaging_time(t_win), wheel_down_s,'k')
            
            hold all;
            pupil_down_s = filter(ones(winsize,1)/winsize, 1, pupil_down);
            pupil_down_s=[pupil_down_s(winsize/2+1:end); zeros(winsize/2,1)];
            
            plot(pupil_time,pupil_Norm/400);
            hold all;
            plot(imaging_time(t_win), pupil_down_s/400,'r')
            legend('wheel','smoothed wheel','pupil','smooth pupil');
            xlabel('Time [samples]');
            hh(2) = subplot(3,1,2);
            
            for kk=1:length(brain_areas)
                %             hh(kk+2)=subplot(length(brain_areas)+2,1,kk+2);
                plot(imaging_time(t_win),TV(ic(kk),:)-(kk-1)*10)
                %         ylabel([ brain_areas{icp(kk)}])
                hold all;
            end
            linkaxes(hh,'x');xlabel('Time [samples]');title('Total Variation');
            axis tight;
            legend(brain_areas,'Location','Best')
            
            
            
            hh(3) = subplot(3,1,3);
            for kk=1:size(cenmat,2)
                %             hh(kk+2)=subplot(length(brain_areas)+2,1,kk+2);
                plot(imaging_time(t_win),squeeze(cenmat(ni,kk,:))-(kk-1)*40)
                %         ylabel([ brain_areas{icp(kk)}])
                hold all;
            end
            linkaxes(hh,'x');xlabel('Time [samples]');title(centralitynames{ni});
            axis tight;
            legend(region_lut(networknums),'Location','Best')
            
        end
        
        
        
        %%5
        [coeff,score,latent] = pca(TV');
        figure;
        plot_by_quentile(TV(ic,:), pupil_down, brain_areas(ic));
        ylim([0 45])
        title('Average Pupil Size');
        figure;
        for ni=1:size(cenmat,1)
            subplot(2,2,ni);
            plot_by_quentile(squeeze(cenmat(ni, networknums,:)), pupil_down, region_lut(networknums));
            title(centralitynames{ni});ylim([0 45])
        end
    end
end

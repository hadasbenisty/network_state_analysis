
function [Trials,TrialsFull]=getTrialAllEvents_wheelonly(outputPath,colorname,dFoF,dFoF_parcells,timestamps,fsimaging,fsspike2,channels_data,Idx)

%% locomotion onset
if ~isempty(timestamps.wheelOn)
    disp(['Extracting locomotion trials in wavelength ',colorname]);
    timestamps.wheelOn=timestamps.wheelOn(timestamps.wheelOn>timestamps.timaging(1)+preEventWin & timestamps.wheelOn<timestamps.timaging(end)+postEventWin); %only get wheel on events during imaging
    
    if size(dFoF_parcells.parcels_time_trace,2)==length(timestamps.timaging)
    elseif size(dFoF_parcells.parcels_time_trace,2)<length(timestamps.timaging)
        timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.parcels_time_trace,2));
    elseif size(dFoF_parcells.parcels_time_trace,2)>length(timestamps.timaging)
        dFoF_parcells.parcels_time_trace=dFoF_parcells.parcels_time_trace(:,1:length(timestamps.timaging));
    end
    
    %removed the line below 08/08/2020
    %timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.parcels_time_trace,2));%remove excess timestamps if you don't have corresponding frames
    [Trials.wheelOn,TrialsFull.wheelOn]=getEventTrials(timestamps.wheelOn,preEventWin,postEventWin,baselineWin,dFoF_parcells.parcels_time_trace,[],timestamps.timaging,fsimaging,fsspike2,channels_data.wheelspeed,Idx,[],[],[],[]);
    tstamps=nanmean(TrialsFull.wheelOn.fulltimeStamps,2);
    % h3=figure;subplot(2,1,1),plot(tstamps,nanmean(Trials.wheelOn.Norm.RightFrontal_dffIntp,2));
    % hold on;plot(tstamps,nanmean(Trials.wheelOn.Norm.RightVis_dffIntp,2));
    % plot(tstamps,nanmean(Trials.wheelOn.Norm.RightSomat_dffIntp,2));
    % plot(tstamps,nanmean(Trials.wheelOn.Norm.LeftFrontal_dffIntp,2));
    % plot(tstamps,nanmean(Trials.wheelOn.Norm.LeftVis_dffIntp,2));
    % plot(tstamps,nanmean(Trials.wheelOn.Norm.LeftSomat_dffIntp,2));title(strcat('WheelOn-Norm-',colorname));legend('R-Frontal','R-Vis','R-Somat','L-Frontal','L-Vis','L-Somat');hold off
    % subplot(2,1,2); plot(-preEventWin:1/fsspike2:(postEventWin-1/fsspike2),nanmean(Trials.wheelOn.RunSpeed,2));ylabel('Speed(cm/s'), xlabel('Time');
    % saveas(h3,fullfile(outputPath,strcat('WheelOn-Norm',colorname)));
    
    
    h4=figure;subplot(2,1,1),plot(tstamps,nanmean(Trials.wheelOn.NonNorm.RightFrontal_dffIntp,2));ylim([-0.02 0.01]);
    hold on;plot(tstamps,nanmean(Trials.wheelOn.NonNorm.RightVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.wheelOn.NonNorm.RightSomat_dffIntp,2));
    plot(tstamps,nanmean(Trials.wheelOn.NonNorm.LeftFrontal_dffIntp,2));
    plot(tstamps,nanmean(Trials.wheelOn.NonNorm.LeftVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.wheelOn.NonNorm.LeftSomat_dffIntp,2));title(strcat('WheelOn-NonNorm-',colorname));ylabel('DFF');legend('R-Frontal','R-Vis','R-Somat','L-Frontal','L-Vis','L-Somat');hold off
    subplot(2,1,2); plot(-plotpreEventWin:1/fsspike2:(postEventWin-1/fsspike2),nanmean(Trials.wheelOn.RunSpeed,2)*100);ylabel('Speed(cm/s'), ylim([-2 20]),xlabel('Time');
    saveas(h4,fullfile(outputPath,strcat('WheelOn-NonNorm_v2',colorname)));
    hold off;
    
    h5=figure;subplot(2,1,1),plot(tstamps,nanmean(Trials.wheelOn.ZNorm.RightFrontal_dffIntp,2));ylim([-2 4]);
    hold on;plot(tstamps,nanmean(Trials.wheelOn.ZNorm.RightVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.wheelOn.ZNorm.RightSomat_dffIntp,2));
    plot(tstamps,nanmean(Trials.wheelOn.ZNorm.LeftFrontal_dffIntp,2));
    plot(tstamps,nanmean(Trials.wheelOn.ZNorm.LeftVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.wheelOn.ZNorm.LeftSomat_dffIntp,2));title(strcat('WheelOn-ZNorm-',colorname)); ylabel('Z-Value');legend('R-Frontal','R-Vis','R-Somat','L-Frontal','L-Vis','L-Somat');hold off
    subplot(2,1,2); plot(-plotpreEventWin:1/fsspike2:(postEventWin-1/fsspike2),nanmean(Trials.wheelOn.RunSpeed,2)*100);ylabel('Speed(cm/s'), ylim([-2 20]),xlabel('Time');
    saveas(h5,fullfile(outputPath,strcat('WheelOn-ZNorm_v2',colorname)));
    hold off;
    % make movie
%     
%     filename=strcat('WheelOn-ZNorm_v2',colorname);
%     videoName=fullfile(outputPath,filename);
%     video=squeeze(nanmean(TrialsFull.wheelOn.ZNorm.fullDFF(:,:,:),2));
%     video=video';
%     finvideo=reshape(video,R,C,size(video,2));climit=[-0.1 0.5];
%     makeMovie(videoName,finvideo,climit,fsimaging);
%     hold off;
%     
%     % make movie
%     filename=strcat('WheelOn-NonNorm_v2',colorname);
%     videoName=fullfile(outputPath,filename);
%     video=squeeze(nanmean(TrialsFull.wheelOn.NonNorm.fullDFF(:,:,:),2));
%     video=video';
%     finvideo=reshape(video,R,C,size(video,2));climit=[-0.02 0.02];
%     %[finvideo]= smoothMovie(finvideo,R,C,2); %smooth movie
%     makeMovie(videoName,finvideo,climit,fsimaging);
%     hold off;
    
    % pick out key timepoints and make brain colormaps
%    clear Zvideo Zfinvideo
%     Zvideo=squeeze(nanmean(TrialsFull.wheelOn.ZNorm.fullDFF,2));
%     Zvideo=Zvideo';
%     Zfinvideo=reshape(Zvideo,R,C,size(Zvideo,2));
%     timepoint1=(plotpreEventWin-1)*fsimaging; %baseline
%     timepoint2=(plotpreEventWin+0.3)*fsimaging; %1/2 s after stimulus onset
%     timepoint3=(plotpreEventWin+3)*fsimaging;%3s after stim onset
%     cmap = [ 1 1 1 ; parula(200) ];
%     ex1=figure;imagesc(Zfinvideo(:,:,timepoint1));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint1)));
%     ex2=figure;imagesc(Zfinvideo(:,:,timepoint2));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint2)));
%     ex3=figure;imagesc(Zfinvideo(:,:,timepoint3));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint3)));
%     saveas(ex1,fullfile(outputPath,strcat('Example1Time-WheelOn-ZNorm_v2',colorname)));
%     saveas(ex2,fullfile(outputPath,strcat('Example2Time-WheelOn-ZNorm_v2',colorname)));
%     saveas(ex3,fullfile(outputPath,strcat('Example3Time-WheelOn-ZNorm_v2',colorname)));
%     hold off;
    figure;subplot(2,1,1);
    r=sqrt(size(Trials.wheelOn.NonNorm.LeftFrontal_dffIntp,2)-1);
    shadedErrorBar(tstamps,nanmean(Trials.wheelOn.NonNorm.LeftFrontal_dffIntp,2),nanstd(Trials.wheelOn.NonNorm.LeftFrontal_dffIntp,0,2)./r,'lineprops','b');ylim([-0.02 0.01]);
    title(strcat('Frontal wheelOn-NonNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.wheelOn.NonNorm.RightFrontal_dffIntp,2),nanstd(Trials.wheelOn.NonNorm.RightFrontal_dffIntp,0,2)./r,'lineprops','g');xline(0);
    
    hold off;
    subplot(2,1,2);
    shadedErrorBar(tstamps,nanmean(Trials.wheelOn.NonNorm.RightVis_dffIntp,2),nanstd(Trials.wheelOn.NonNorm.RightVis_dffIntp,0,2)./r,'lineprops','b');hold on;ylim([-0.02 0.01]);
    shadedErrorBar(tstamps,nanmean(Trials.wheelOn.NonNorm.LeftVis_dffIntp,2),nanstd(Trials.wheelOn.NonNorm.LeftVis_dffIntp,0,2)./r,'lineprops','g');
    ylim([-0.02 0.01]);title(strcat('Visual wheelOn-NonNorm',colorname));ylabel('dFoF');xlabel('Time(s)');xline(0);
    mysave(gcf, fullfile(outputPath,strcat('Frontal_vs_Vis_wheelOn-NonNorm_v2',colorname)), 'all');
    % filename=strcat('WheelOn-ZScore',colorname);
    % videoName=fullfile(outputPath,filename);
    % video=squeeze(nanmean(TrialsFull.wheelOn.ZNorm.fullDFF,2));
    % video=video';
    % finvideo=reshape(video,R,C,size(video,2));climit=[-2 2];
    % [finvideo]= smoothMovie(finvideo,R,C,2); %smooth movie
    % makeMovie(videoName,finvideo,climit,fsimaging);
    hold off;
    figure;subplot(2,1,1);
    r=sqrt(size(Trials.wheelOn.NonNorm.LeftFrontal_dffIntp,2)-1);
    shadedErrorBar(tstamps,nanmean(Trials.wheelOn.ZNorm.LeftFrontal_dffIntp,2),nanstd(Trials.wheelOn.ZNorm.LeftFrontal_dffIntp,0,2)./r,'lineprops','b');
    title(strcat('Frontal wheelOn-ZNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.wheelOn.ZNorm.RightFrontal_dffIntp,2),nanstd(Trials.wheelOn.ZNorm.RightFrontal_dffIntp,0,2)./r,'lineprops','g');xline(0);
    
    subplot(2,1,2);
    shadedErrorBar(tstamps,nanmean(Trials.wheelOn.ZNorm.RightVis_dffIntp,2),nanstd(Trials.wheelOn.ZNorm.RightVis_dffIntp,0,2)./r,'lineprops','b');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.wheelOn.ZNorm.LeftVis_dffIntp,2),nanstd(Trials.wheelOn.ZNorm.LeftVis_dffIntp,0,2)./r,'lineprops','g');xline(0);
    title(strcat('Visual wheelOn-ZNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    mysave(gcf, fullfile(outputPath,strcat('Frontal_vs_Vis_wheelOn-ZNorm_v2',colorname)), 'all');
    hold off;
end


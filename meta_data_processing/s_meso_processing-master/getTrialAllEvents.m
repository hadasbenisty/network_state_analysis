
function [Trials,TrialsFull]=getTrialAllEvents(params,outputPath,colorname,dFoF,dFoF_parcells,timestamps,fsimaging,fsspike2,channels_data,Idx,R,C)
%% airpuff
plotpreEventWin=2;preEventWin=2; postEventWin=5; baselineWin=2;
if params.airpuffAn&&~isempty(timestamps.airpuff)
    disp(['Extracting airpuff trials in wavelength ',colorname]);
    if size(dFoF_parcells.blue,2)==length(timestamps.timaging)
    elseif size(dFoF_parcells.blue,2)<length(timestamps.timaging)
        timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.blue,2));
    elseif size(dFoF_parcells.blue,2)>length(timestamps.timaging)
        dFoF_parcells.blue=dFoF_parcells.blue(1:length(timestamps.timaging),:);
    end
    %timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.blue,2));
    %removed 08082020
    timestamps.airpuff=timestamps.airpuff(timestamps.airpuff>timestamps.timaging(1)+preEventWin & timestamps.airpuff<timestamps.timaging(end)+postEventWin); %only get airpuff events during imaging
    [Trials.airpuff, TrialsFull.airpuff]=getEventTrials(timestamps.airpuff,preEventWin,postEventWin,baselineWin,dFoF_parcells.(colorname),dFoF.(colorname),timestamps.timaging,fsimaging,fsspike2,channels_data.wheelspeed,Idx,[],[],[],[]);
    tstamps=nanmean(TrialsFull.airpuff.fulltimeStamps,2);
    h=figure;subplot(2,1,1);
    plot(tstamps,nanmean(Trials.airpuff.ZNorm.RightFrontal_dffIntp,2));ylim([-2 4]);
    hold on;plot(tstamps,nanmean(Trials.airpuff.ZNorm.RightVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.airpuff.ZNorm.RightSomat_dffIntp,2));
    plot(tstamps,nanmean(Trials.airpuff.ZNorm.LeftFrontal_dffIntp,2));
    plot(tstamps,nanmean(Trials.airpuff.ZNorm.LeftVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.airpuff.ZNorm.LeftSomat_dffIntp,2));title(strcat('Airpuff-ZNorm',colorname));ylabel('Z-value');xlabel('Time(s)');hold off
    subplot(2,1,2); plot(-plotpreEventWin:1/fsspike2:(postEventWin-1/fsspike2),nanmean(Trials.airpuff.RunSpeed,2)*100);ylabel('Speed(cm/s'), ylim([-2 20]),xlabel('Time');
    saveas(h,fullfile(outputPath,strcat('Airpuff-ZNorm_v2',colorname)))
    hold off;
    
    
    h1=figure;subplot(2,1,1);
    plot(tstamps,nanmean(Trials.airpuff.NonNorm.RightFrontal_dffIntp,2));ylim([-0.02 0.01]);
    hold on;plot(tstamps,nanmean(Trials.airpuff.NonNorm.RightVis_dffIntp,2));ylim([-0.02 0.01]);
    plot(tstamps,nanmean(Trials.airpuff.NonNorm.RightSomat_dffIntp,2));ylim([-0.02 0.01]);
    plot(tstamps,nanmean(Trials.airpuff.NonNorm.LeftFrontal_dffIntp,2));ylim([-0.02 0.01]);
    plot(tstamps,nanmean(Trials.airpuff.NonNorm.LeftVis_dffIntp,2));ylim([-0.02 0.01]);
    plot(tstamps,nanmean(Trials.airpuff.NonNorm.LeftSomat_dffIntp,2));title(strcat('Airpuff-NonNorm',colorname));ylabel('dff');xlabel('Time(s)');hold off
    subplot(2,1,2); plot(-plotpreEventWin:1/fsspike2:(postEventWin-1/fsspike2),nanmean(Trials.airpuff.RunSpeed,2)*100);ylabel('Speed(cm/s'), ylim([-2 20]),xlabel('Time');
    saveas(h1,fullfile(outputPath,strcat('Airpuff-NonNorm_v2',colorname)));legend('R-Frontal','R-Vis','R-Somat','L-Frontal','L-Vis','L-Somat');
    %legend('R-Frontal','R-Vis','R-Somat','L-Frontal','L-Vis','L-Somat');
    hold off;
    
    figure;subplot(2,1,1);
    r=sqrt(size(Trials.airpuff.NonNorm.LeftFrontal_dffIntp,2)-1);
    shadedErrorBar(tstamps,nanmean(Trials.airpuff.NonNorm.LeftFrontal_dffIntp,2),nanstd(Trials.airpuff.NonNorm.LeftFrontal_dffIntp,0,2)./r,'lineprops','b');ylim([-0.02 0.01]);
    title(strcat('Frontal Airpuff-NonNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.airpuff.NonNorm.RightFrontal_dffIntp,2),nanstd(Trials.airpuff.NonNorm.RightFrontal_dffIntp,0,2)./r,'lineprops','g');xline(0);
    
    subplot(2,1,2);
    shadedErrorBar(tstamps,nanmean(Trials.airpuff.NonNorm.RightVis_dffIntp,2),nanstd(Trials.airpuff.NonNorm.RightVis_dffIntp,0,2)./r,'lineprops','b');hold on;ylim([-0.02 0.01])
    shadedErrorBar(tstamps,nanmean(Trials.airpuff.NonNorm.LeftVis_dffIntp,2),nanstd(Trials.airpuff.NonNorm.LeftVis_dffIntp,0,2)./r,'lineprops','g');xline(0);
    ylim([-0.02 0.01]);title(strcat('Visual Airpuff-NonNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    mysave(gcf, fullfile(outputPath,strcat('Frontal_vs_Vis_Airpuff-NonNorm_v2',colorname)), 'all');
    hold off;
    
    figure;subplot(2,1,1);
    r=sqrt(size(Trials.airpuff.NonNorm.LeftFrontal_dffIntp,2)-1);
    shadedErrorBar(tstamps,nanmean(Trials.airpuff.ZNorm.LeftFrontal_dffIntp,2),nanstd(Trials.airpuff.ZNorm.LeftFrontal_dffIntp,0,2)./r,'lineprops','b');
    title(strcat('Frontal Airpuff-ZNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.airpuff.ZNorm.RightFrontal_dffIntp,2),nanstd(Trials.airpuff.ZNorm.RightFrontal_dffIntp,0,2)./r,'lineprops','g');xline(0);
    
    subplot(2,1,2);
    shadedErrorBar(tstamps,nanmean(Trials.airpuff.ZNorm.RightVis_dffIntp,2),nanstd(Trials.airpuff.ZNorm.RightVis_dffIntp,0,2)./r,'lineprops','g');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.airpuff.ZNorm.LeftVis_dffIntp,2),nanstd(Trials.airpuff.ZNorm.LeftVis_dffIntp,0,2)./r,'lineprops','b');xline(0);
    title(strcat('Visual Airpuff-ZNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    mysave(gcf, fullfile(outputPath,strcat('Frontal_vs_Vis_Airpuff-ZNorm_v2',colorname)), 'all');
    hold off;
    
    % make movie
    % filename=strcat('Airpuff-Norm',colorname);
    % videoName=fullfile(outputPath,filename);
    % video=squeeze(nanmean(TrialsFull.airpuff.Norm.fullDFF,2));
    % video=video';
    % finvideo=reshape(video,R,C,size(video,2));climit=[-0.1 0.3];
    % makeMovie(videoName,finvideo,climit,fsimaging);
    hold off;
    filename=strcat('Airpuff-NonNorm',colorname);
    videoName=fullfile(outputPath,filename);
    video=squeeze(nanmean(TrialsFull.airpuff.NonNorm.fullDFF,2));
    video=video';
    finvideo=reshape(video,R,C,size(video,2));climit=[-0.02 0.02];
    %[finvideo]= smoothMovie(finvideo,R,C,2); %smooth movie
    makeMovie(videoName,finvideo,climit,fsimaging);
    
    % pick out key timepoints and make brain colormaps
    Zvideo=squeeze(nanmean(TrialsFull.airpuff.ZNorm.fullDFF,2));
    Zvideo=Zvideo';
    Zfinvideo=reshape(Zvideo,R,C,size(Zvideo,2));
    timepoint1=(plotpreEventWin-1)*fsimaging; %baseline
    timepoint2=(plotpreEventWin+0.3)*fsimaging; %1/2 s after stimulus onset
    timepoint3=(plotpreEventWin+3)*fsimaging;%3s after stim onset
    cmap = [ 1 1 1 ; parula(200) ];
    ex1=figure;x=imagesc(Zfinvideo(:,:,timepoint1)); colormap(cmap);caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint1)));
    ex2=figure;imagesc(Zfinvideo(:,:,timepoint2));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint2)));
    ex3=figure;imagesc(Zfinvideo(:,:,timepoint3));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint3)));
    saveas(ex1,fullfile(outputPath,strcat('Example1Time-Airpuff-ZNorm_v2',colorname)));
    saveas(ex2,fullfile(outputPath,strcat('Example2Time-Airpuff-ZNorm_v2',colorname)));
    saveas(ex3,fullfile(outputPath,strcat('Example3Time-Airpuff-ZNorm_v2',colorname)));
else
end

%% locomotion onset
if ~isempty(timestamps.wheelOn)
    disp(['Extracting locomotion trials in wavelength ',colorname]);
    timestamps.wheelOn=timestamps.wheelOn(timestamps.wheelOn>timestamps.timaging(1)+preEventWin & timestamps.wheelOn<timestamps.timaging(end)+postEventWin); %only get wheel on events during imaging
    
    if size(dFoF_parcells.blue,2)==length(timestamps.timaging)
    elseif size(dFoF_parcells.blue,2)<length(timestamps.timaging)
        timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.blue,2));
    elseif size(dFoF_parcells.blue,2)>length(timestamps.timaging)
        dFoF_parcells.blue=dFoF_parcells.blue(:,1:length(timestamps.timaging));
    end
    
    %removed the line below 08/08/2020
    %timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.blue,2));%remove excess timestamps if you don't have corresponding frames
    [Trials.wheelOn,TrialsFull.wheelOn]=getEventTrials(timestamps.wheelOn,preEventWin,postEventWin,baselineWin,dFoF_parcells.(colorname),dFoF.(colorname),timestamps.timaging,fsimaging,fsspike2,channels_data.wheelspeed,Idx,[],[],[],[]);
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
    
    filename=strcat('WheelOn-ZNorm_v2',colorname);
    videoName=fullfile(outputPath,filename);
    video=squeeze(nanmean(TrialsFull.wheelOn.ZNorm.fullDFF(:,:,:),2));
    video=video';
    finvideo=reshape(video,R,C,size(video,2));climit=[-0.1 0.5];
    makeMovie(videoName,finvideo,climit,fsimaging);
    hold off;
    
    % make movie
    filename=strcat('WheelOn-NonNorm_v2',colorname);
    videoName=fullfile(outputPath,filename);
    video=squeeze(nanmean(TrialsFull.wheelOn.NonNorm.fullDFF(:,:,:),2));
    video=video';
    finvideo=reshape(video,R,C,size(video,2));climit=[-0.02 0.02];
    %[finvideo]= smoothMovie(finvideo,R,C,2); %smooth movie
    makeMovie(videoName,finvideo,climit,fsimaging);
    hold off;
    
    % pick out key timepoints and make brain colormaps
    clear Zvideo Zfinvideo
    Zvideo=squeeze(nanmean(TrialsFull.wheelOn.ZNorm.fullDFF,2));
    Zvideo=Zvideo';
    Zfinvideo=reshape(Zvideo,R,C,size(Zvideo,2));
    timepoint1=(plotpreEventWin-1)*fsimaging; %baseline
    timepoint2=(plotpreEventWin+0.3)*fsimaging; %1/2 s after stimulus onset
    timepoint3=(plotpreEventWin+3)*fsimaging;%3s after stim onset
    cmap = [ 1 1 1 ; parula(200) ];
    ex1=figure;imagesc(Zfinvideo(:,:,timepoint1));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint1)));
    ex2=figure;imagesc(Zfinvideo(:,:,timepoint2));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint2)));
    ex3=figure;imagesc(Zfinvideo(:,:,timepoint3));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint3)));
    saveas(ex1,fullfile(outputPath,strcat('Example1Time-WheelOn-ZNorm_v2',colorname)));
    saveas(ex2,fullfile(outputPath,strcat('Example2Time-WheelOn-ZNorm_v2',colorname)));
    saveas(ex3,fullfile(outputPath,strcat('Example3Time-WheelOn-ZNorm_v2',colorname)));
    hold off;
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
    
    if params.visStimAn && ~isempty(timestamps.visstart)
        
    disp(['Extracting visual trials in wavelength ',colorname]);
    timestamps.visstart=timestamps.visstart(timestamps.visstart>timestamps.timaging(1)+preEventWin & timestamps.visstart<timestamps.timaging(end)+postEventWin); %only get wheel on events during imaging
    
    timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.blue,2));%remove excess timestamps if you don't have corresponding frames
    [Trials.visstart,TrialsFull.visstart]=getEventTrials(timestamps.visstart,preEventWin,postEventWin,baselineWin,dFoF_parcells.(colorname),dFoF.(colorname),timestamps.timaging,fsimaging,fsspike2,channels_data.wheelspeed,Idx,[],[],[],[]);
    tstamps=nanmean(TrialsFull.visstart.fulltimeStamps,2);
    % h3=figure;subplot(2,1,1),plot(tstamps,nanmean(Trials.visstart.Norm.RightFrontal_dffIntp,2));
    % hold on;plot(tstamps,nanmean(Trials.visstart.Norm.RightVis_dffIntp,2));
    % plot(tstamps,nanmean(Trials.visstart.Norm.RightSomat_dffIntp,2));
    % plot(tstamps,nanmean(Trials.visstart.Norm.LeftFrontal_dffIntp,2));
    % plot(tstamps,nanmean(Trials.visstart.Norm.LeftVis_dffIntp,2));
    % plot(tstamps,nanmean(Trials.visstart.Norm.LeftSomat_dffIntp,2));title(strcat('visstart-Norm-',colorname));legend('R-Frontal','R-Vis','R-Somat','L-Frontal','L-Vis','L-Somat');hold off
    % subplot(2,1,2); plot(-preEventWin:1/fsspike2:(postEventWin-1/fsspike2),nanmean(Trials.visstart.RunSpeed,2));ylabel('Speed(cm/s'), xlabel('Time');
    % saveas(h3,fullfile(outputPath,strcat('visstart-Norm',colorname)));
    
    
    h4=figure;subplot(2,1,1),plot(tstamps,nanmean(Trials.visstart.NonNorm.RightFrontal_dffIntp,2));;ylim([-0.02 0.01]);
    hold on;plot(tstamps,nanmean(Trials.visstart.NonNorm.RightVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.visstart.NonNorm.RightSomat_dffIntp,2));
    plot(tstamps,nanmean(Trials.visstart.NonNorm.LeftFrontal_dffIntp,2));
    plot(tstamps,nanmean(Trials.visstart.NonNorm.LeftVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.visstart.NonNorm.LeftSomat_dffIntp,2));title(strcat('visstart-NonNorm-',colorname));ylabel('DFF');legend('R-Frontal','R-Vis','R-Somat','L-Frontal','L-Vis','L-Somat');hold off
    subplot(2,1,2); plot(-plotpreEventWin:1/fsspike2:(postEventWin-1/fsspike2),nanmean(Trials.visstart.RunSpeed,2)*100);ylabel('Speed(cm/s'), ylim([-2 20]),xlabel('Time');
    saveas(h4,fullfile(outputPath,strcat('visstart-NonNorm_v2',colorname)));
    hold off;
    
    h5=figure;subplot(2,1,1),plot(tstamps,nanmean(Trials.visstart.ZNorm.RightFrontal_dffIntp,2));ylim([-2 4]);
    hold on;plot(tstamps,nanmean(Trials.visstart.ZNorm.RightVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.visstart.ZNorm.RightSomat_dffIntp,2));
    plot(tstamps,nanmean(Trials.visstart.ZNorm.LeftFrontal_dffIntp,2));
    plot(tstamps,nanmean(Trials.visstart.ZNorm.LeftVis_dffIntp,2));
    plot(tstamps,nanmean(Trials.visstart.ZNorm.LeftSomat_dffIntp,2));title(strcat('visstart-ZNorm-',colorname)); ylabel('Z-Value');legend('R-Frontal','R-Vis','R-Somat','L-Frontal','L-Vis','L-Somat');hold off
    subplot(2,1,2); plot(-plotpreEventWin:1/fsspike2:(postEventWin-1/fsspike2),nanmean(Trials.visstart.RunSpeed,2)*100);ylabel('Speed(cm/s'), ylim([-2 20]),xlabel('Time');
    saveas(h5,fullfile(outputPath,strcat('visstart-ZNorm_v2',colorname)));
    hold off;
    % make movie
    
    filename=strcat('visstart-ZNorm_v2',colorname);
    videoName=fullfile(outputPath,filename);
    video=squeeze(nanmean(TrialsFull.visstart.ZNorm.fullDFF(:,:,:),2));
    video=video';
    finvideo=reshape(video,R,C,size(video,2));climit=[-0.1 0.5];
    makeMovie(videoName,finvideo,climit,fsimaging);
    hold off;
    
    % make movie
    filename=strcat('visstart-NonNorm_v2',colorname);
    videoName=fullfile(outputPath,filename);
    video=squeeze(nanmean(TrialsFull.visstart.NonNorm.fullDFF(:,:,:),2));
    video=video';
    finvideo=reshape(video,R,C,size(video,2));climit=[-0.02 0.02];
    %[finvideo]= smoothMovie(finvideo,R,C,2); %smooth movie
    makeMovie(videoName,finvideo,climit,fsimaging);
    hold off;
    
    % pick out key timepoints and make brain colormaps
    clear Zvideo Zfinvideo
    Zvideo=squeeze(nanmean(TrialsFull.visstart.ZNorm.fullDFF,2));
    Zvideo=Zvideo';
    Zfinvideo=reshape(Zvideo,R,C,size(Zvideo,2));
    timepoint1=(plotpreEventWin-1)*fsimaging; %baseline
    timepoint2=(plotpreEventWin+0.3)*fsimaging; %1/2 s after stimulus onset
    timepoint3=(plotpreEventWin+3)*fsimaging;%3s after stim onset
    cmap = [ 1 1 1 ; parula(200) ];
    ex1=figure;imagesc(Zfinvideo(:,:,timepoint1));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint1)));
    ex2=figure;imagesc(Zfinvideo(:,:,timepoint2));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint2)));
    ex3=figure;imagesc(Zfinvideo(:,:,timepoint3));colormap(cmap); caxis([-1 5]);colorbar;title(strcat('Timepoint',num2str(timepoint3)));
    saveas(ex1,fullfile(outputPath,strcat('Example1Time-visstart-ZNorm_v2',colorname)));
    saveas(ex2,fullfile(outputPath,strcat('Example2Time-visstart-ZNorm_v2',colorname)));
    saveas(ex3,fullfile(outputPath,strcat('Example3Time-visstart-ZNorm_v2',colorname)));
    hold off;
    figure;subplot(2,1,1);
    r=sqrt(size(Trials.visstart.NonNorm.LeftFrontal_dffIntp,2)-1);
    shadedErrorBar(tstamps,nanmean(Trials.visstart.NonNorm.LeftFrontal_dffIntp,2),nanstd(Trials.visstart.NonNorm.LeftFrontal_dffIntp,0,2)./r,'lineprops','b');ylim([-0.02 0.01]);
    title(strcat('Frontal visstart-NonNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.visstart.NonNorm.RightFrontal_dffIntp,2),nanstd(Trials.visstart.NonNorm.RightFrontal_dffIntp,0,2)./r,'lineprops','g');xline(0);
    
    hold off;
    subplot(2,1,2);
    shadedErrorBar(tstamps,nanmean(Trials.visstart.NonNorm.RightVis_dffIntp,2),nanstd(Trials.visstart.NonNorm.RightVis_dffIntp,0,2)./r,'lineprops','b');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.visstart.NonNorm.LeftVis_dffIntp,2),nanstd(Trials.visstart.NonNorm.LeftVis_dffIntp,0,2)./r,'lineprops','g');
    ylim([-0.02 0.01]);title(strcat('Visual visstart-NonNorm',colorname));ylabel('dFoF');xlabel('Time(s)');xline(0);
    mysave(gcf, fullfile(outputPath,strcat('Frontal_vs_Vis_visstart-NonNorm_v2',colorname)), 'all');hold off;
    

    figure;subplot(2,1,1);
    r=sqrt(size(Trials.visstart.NonNorm.LeftFrontal_dffIntp,2)-1);
    shadedErrorBar(tstamps,nanmean(Trials.visstart.ZNorm.LeftFrontal_dffIntp,2),nanstd(Trials.visstart.ZNorm.LeftFrontal_dffIntp,0,2)./r,'lineprops','b');ylim([-2 4]);
    title(strcat('Frontal visstart-ZNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    shadedErrorBar(tstamps,nanmean(Trials.visstart.ZNorm.RightFrontal_dffIntp,2),nanstd(Trials.visstart.ZNorm.RightFrontal_dffIntp,0,2)./r,'lineprops','g');xline(0);
    
    subplot(2,1,2);
    shadedErrorBar(tstamps,nanmean(Trials.visstart.ZNorm.RightVis_dffIntp,2),nanstd(Trials.visstart.ZNorm.RightVis_dffIntp,0,2)./r,'lineprops','g');hold on;ylim([-2 4]);
    shadedErrorBar(tstamps,nanmean(Trials.visstart.ZNorm.LeftVis_dffIntp,2),nanstd(Trials.visstart.ZNorm.LeftVis_dffIntp,0,2)./r,'lineprops','b');xline(0);
    title(strcat('Visual visstart-ZNorm',colorname));ylabel('dFoF');xlabel('Time(s)');hold on;
    mysave(gcf, fullfile(outputPath,strcat('Frontal_vs_Vis_visstart-ZNorm_v2',colorname)), 'all');
    hold off;
    else
    end
end


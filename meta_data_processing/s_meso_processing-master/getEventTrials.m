function [Trial, TrialFull]=getEventTrials(eventTS,preEventWin,postEventWin,baselineWin,final_dff_parcellated,final_dff,imagingTS,fsimaging,fsspike2,RunSpeed,Idx,EEG,proc,pupilCamTimes,fspupilcam)
%if ~isempty(eventTS)
    eventTS=eventTS(:)';
    [Trial.NonNorm.parcellsDFF,Trial.parcellstimeStamps]= TrialTimeArrangeDff(final_dff_parcellated,imagingTS,fsimaging,preEventWin,eventTS,postEventWin);%do parcells only
    %Trial.NonNorm.parcellsDFF=fillmissing(Trial.NonNorm.parcellsDFF,'linear',1,'EndValues','nearest');%interpolate
    %line below commented out jan 7 2021 lav
    %[TrialFull.NonNorm.fullDFF,TrialFull.fulltimeStamps]= TrialTimeArrangeDff(final_dff,imagingTS,fsimaging,preEventWin,eventTS,postEventWin);% do the full movie
    %Trial.NonNorm.fullDFF=fillmissing(Trial.NonNorm.fullDFF,'linear',1,'EndValues','nearest');%interpolate
    
    %normalize by baseline before event,parcellated
    firstTimeStamp=Trial.parcellstimeStamps(1,1);
    BaselineIdx=Trial.parcellstimeStamps(:,1)>=firstTimeStamp & Trial.parcellstimeStamps(:,1)<(firstTimeStamp+preEventWin);
    Baseline_Parcells=Trial.NonNorm.parcellsDFF(BaselineIdx,:,:);
    MeanBaseline_P=nanmean(Baseline_Parcells,1); %normalize with each trialmean
    StdBaseline_P=nanstd(Baseline_Parcells,1);
    %MeanBaseline_P=nanmean(nanmean(Baseline_Parcells,1)); %normalize with mean baseline across trials
    %StdBaseline_P=nanstd(MeanBaseline_P,1); %std of mean baseline across triasl
    MeanSubBaseline_P=repmat(MeanBaseline_P,size(Trial.NonNorm.parcellsDFF,1),1,1);
    tmp1=bsxfun(@minus,Trial.NonNorm.parcellsDFF,MeanSubBaseline_P);
    tmp2=bsxfun(@rdivide,tmp1,MeanSubBaseline_P);
    Trial.Norm.parcellsDFF=tmp2; %dff with f0 as pres-stimulus period for each trial
    Trial.ZNorm.parcellsDFF=(Trial.NonNorm.parcellsDFF-MeanBaseline_P)./StdBaseline_P;
    
    %normalize by baseline before event,raw non-parcellated
    firstTimeStamp=TrialFull.fulltimeStamps(1,1);
    BaselineIdx=TrialFull.fulltimeStamps(:,1)>=firstTimeStamp & TrialFull.fulltimeStamps(:,1)<(firstTimeStamp+baselineWin);
    %BaselineIdx=Trial.fulltimeStamps(:,1)<0;
    Baseline_Full=TrialFull.NonNorm.fullDFF(BaselineIdx,:,:);
    MeanBaseline=nanmean(Baseline_Full,1);
    %MeanBaseline=nanmean(nanmean(Baseline_Full,1)); %normalize with mean baseline across trials
    StdBaseline=nanstd(Baseline_Full,1);
    %StdBaseline=nanstd(MeanBaseline,1); %std of mean baseline across triasl
    MeanSubBaseline=repmat(MeanBaseline,size(TrialFull.NonNorm.fullDFF,1),1,1);
    tmp3=bsxfun(@minus,TrialFull.NonNorm.fullDFF,MeanSubBaseline);
    tmp4=bsxfun(@rdivide,tmp3,MeanSubBaseline);
    TrialFull.Norm.fullDFF=tmp4; %dff with f0 as pres-stimulus period for each
    TrialFull.ZNorm.fullDFF=(TrialFull.NonNorm.fullDFF-MeanBaseline)./StdBaseline;
    
    % combined across parcells, non normalized
    Trial.NonNorm.LeftVis_dffIntp=nanmean(Trial.NonNorm.parcellsDFF(:,:,Idx.Left_Vis_idx),3);
    Trial.NonNorm.RightVis_dffIntp=nanmean(Trial.NonNorm.parcellsDFF(:,:,Idx.Right_Vis_idx),3);
    Trial.NonNorm.LeftSomat_dffIntp=nanmean(Trial.NonNorm.parcellsDFF(:,:,Idx.Left_Somat_idx),3);
    Trial.NonNorm.RightSomat_dffIntp=nanmean(Trial.NonNorm.parcellsDFF(:,:,Idx.Right_Somat_idx),3);
    Trial.NonNorm.LeftFrontal_dffIntp=nanmean(Trial.NonNorm.parcellsDFF(:,:,Idx.Left_Frontal_idx),3);
    Trial.NonNorm.RightFrontal_dffIntp=nanmean(Trial.NonNorm.parcellsDFF(:,:,Idx.Right_Frontal_idx),3);
    
    % combined across parcells, normalized  to pre-event baseline, df/f
    Trial.Norm.LeftVis_dffIntp=nanmean(Trial.Norm.parcellsDFF(:,:,Idx.Left_Vis_idx),3);
    Trial.Norm.RightVis_dffIntp=nanmean(Trial.Norm.parcellsDFF(:,:,Idx.Right_Vis_idx),3);
    Trial.Norm.LeftSomat_dffIntp=nanmean(Trial.Norm.parcellsDFF(:,:,Idx.Left_Somat_idx),3);
    Trial.Norm.RightSomat_dffIntp=nanmean(Trial.Norm.parcellsDFF(:,:,Idx.Right_Somat_idx),3);
    Trial.Norm.LeftFrontal_dffIntp=nanmean(Trial.Norm.parcellsDFF(:,:,Idx.Left_Frontal_idx),3);
    Trial.Norm.RightFrontal_dffIntp=nanmean(Trial.Norm.parcellsDFF(:,:,Idx.Right_Frontal_idx),3);
    
    % combined across parcells, ZNormalized  to pre-event baseline, z-score
    Trial.ZNorm.LeftVis_dffIntp=nanmean(Trial.ZNorm.parcellsDFF(:,:,Idx.Left_Vis_idx),3);
    Trial.ZNorm.RightVis_dffIntp=nanmean(Trial.ZNorm.parcellsDFF(:,:,Idx.Right_Vis_idx),3);
    Trial.ZNorm.LeftSomat_dffIntp=nanmean(Trial.ZNorm.parcellsDFF(:,:,Idx.Left_Somat_idx),3);
    Trial.ZNorm.RightSomat_dffIntp=nanmean(Trial.ZNorm.parcellsDFF(:,:,Idx.Right_Somat_idx),3);
    Trial.ZNorm.LeftFrontal_dffIntp=nanmean(Trial.ZNorm.parcellsDFF(:,:,Idx.Left_Frontal_idx),3);
    Trial.ZNorm.RightFrontal_dffIntp=nanmean(Trial.ZNorm.parcellsDFF(:,:,Idx.Right_Frontal_idx),3);
    
    %get run speed
    Trial.RunSpeed=nan((preEventWin+postEventWin)*fsspike2,size(Trial.parcellstimeStamps,2));
    spike2Time=1/fsspike2:1/fsspike2:(length(RunSpeed)/fsspike2);
    
    for i=1:size(Trial.parcellstimeStamps,2)
        tmp=RunSpeed(find(spike2Time>=(eventTS(i)-preEventWin) & spike2Time<=(eventTS(i)+postEventWin)));
        Trial.RunSpeed(:,i)=tmp(1:size(Trial.RunSpeed,1));
    end
%else
%end
end

% %extract pupil and face data for each trial
% if exist('proc','var') && ~isempty('proc')
%     %get pupil
%     pupil=proc.pupil.area';
%     [Trial.pupil,Trial.pupiltimeStamps]= TrialTimeArrangeDff(pupil,pupilCamTimes,fspupilcam,preEventWin,eventTS,postEventWin);
%     if isfield(proc, 'motSVD')
%         %get face PCs
%         wholeFaceSVD=proc.motSVD{1,1};
%         for i=1:size(wholeFaceSVD,2)
%             currpca=wholeFaceSVD(:,i)';
%             [Trial.wholeFaceSVD{i},Trial.wholeFacetimeStamps{i}]= TrialTimeArrangeDff(currpca,pupilCamTimes,fspupilcam,preEventWin,eventTS,postEventWin);
%
%         end
%
%         %get whisker PCs
%         whiskerSVD=proc.motSVD{1,2};
%         for i=1:size(whiskerSVD,2)
%             currpca=whiskerSVD(:,i)';
%             [Trial.whiskerSVD{i},Trial.whiskertimeStamps{i}]= TrialTimeArrangeDff(currpca,pupilCamTimes,fspupilcam,preEventWin,eventTS,postEventWin);
%         end
%     end
% end



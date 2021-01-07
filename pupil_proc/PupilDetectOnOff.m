function[PupOnTS, PupOffTS]= PupilDetectOnOff(data,threshold1,threshold2,stimITI,spike2SR)

visTrace = medfilt1(data,100); %
vis=tsmovavg(visTrace,'s',50,1); %
vis=vis-nanmean(vis);
vis=(vis/max(vis)>threshold1); %use 0.6 threshold to get ride of start and end artifacts
[vison,visoff]=squaredetect(vis,threshold2);

PupOnTS=vison*1/spike2SR;
PupOffTS=visoff*1/spike2SR;
tmp4=diff(PupOnTS); 
tmp5=find(tmp4>(stimITI+0.01)| tmp4<(stimITI-0.01)); %remove artifact in on /off times if the iti time is greater or less than the said iti
warning('Way too many pupil triggers detected, correcting  but double check the output'); 
if ~isempty(tmp5)
    PupOnTS(tmp5)=[];
    PupOffTS(tmp5)=[];
end
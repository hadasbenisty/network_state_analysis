% Load in all the variables for each single mouse from step1, (meso frameindices for each day)
% align baselined and Z-scored Ca/ running/ pupil/ eyeblink traces by running onsets/ running offsets (2+2s)/ spontaneous
% blinks (2+2s)/ whisking onset/ visual cue onsets (4+5.5s) / random timepoints (9.5s)
% save in '_caaligned.mat' in each animal's folder
clear all;
clc;
addpath(genpath('utils'));

% identify data to analyze
% load in log file col1-animalid; col2-task (Cond or Disc); 3-Target (ST or SC or PN or NS);
% col4- Layer(L2 or L5); col5-Useful; col6-RF mapping days; col7-CRF/ORF test before
% learning days; col8-Airpuff CS+ CS- unpaired days; col9-Learning days; col10-Learning plateau;
% col11-CS+ shift days; col12-CS+ shift days plateau; col13- Psychtest days; col14- Psychtestc50 days;
% col15- CRF/ORF after learning days; col16- extinction days; col17- Airpuff CS+ CS- unpaired days after training;
% col18- CS+ orientation; col19- CS- orientation; col20~22-Learning Phase1~3;

%pathname=strcat('\\128.36.220.167\vivo\Lan\Meso-imaging');
pathname=strcat('\\128.36.220.173\vivo2\Lan\Meso-imaging');

[~,animallog]=xlsread(strcat(pathname,'\animal_log_mesoimaging_lt.xlsx'));
% now determine what animal and what column to analyze
target='Meso';
layer='s';
load('\\128.36.220.173\vivo2\Lan\Meso-imaging\template\subregion.mat','allregionspix','allregions');
npix=length(allregionspix); % pixel numbers defining the brain regionl
analysisType = 8;
[fieldname, savesuffix, col] = get_preprocessed_files_suffix(analysisType);


%animalnum=find(strcmp(animallog(:,5),'N')~=1 & strcmp(animallog(:,3),target)==1 & strcmp(animallog(:,4),layer)==1 & strcmp(animallog(:,col),'[]')~=1);
%animalnum=find(strcmp(animallog(:,2),'Cond')==1 & strcmp(animallog(:,5),'N')~=1 & strcmp(animallog(:,col),'[]')~=1);
animalnum=find(strcmp(animallog(:,2),'Cond')==1 & strcmp(animallog(:,5),'N')~=1 & strcmp(animallog(:,col),'[]')~=1);
nanimal=length(animalnum);
animalid=animallog(animalnum,1);

analwin=[];
analwin.win_sti_before=-2.5;
analwin.win_sti_after=3;
analwin.win_loco_before=-2.5;
analwin.win_loco_after=3;

blwin=[];
blwin.win_sti_st=-1;
blwin.win_sti_ed=0;
blwin.win_loco_st=-2;
blwin.win_loco_ed=-1;

if ~isempty(animalnum)
    for rowind=1:nanimal
        animali=animalnum(rowind);
        
        savefilename=fullfile(pathname,char(animalid(rowind)),strcat(char(animalid(rowind)),savesuffix));
        
        display(strcat('loading in file:',savefilename));
        load(strcat(savefilename,'.mat'));
        caframeoffset=floor(caframerate/10);
        % split caframeind by day
        caframeind_split=[0;find(diff(caframeind_all)<0);length(caframeind_all)];
        caframeind_day=eval(char(animallog(animali,col)));
        % set up cropping windows
        spike2rate=spike2rate./downsamplefac;
        %
        winframe=[];
        winframe.sti_st=ceil(-analwin.win_sti_before*caframerate);
        winframe.loco_st=ceil(-analwin.win_loco_before*caframerate);
        winframe.blwin_sti_st=max(ceil((-analwin.win_sti_before+blwin.win_sti_st)*caframerate),1);
        winframe.blwin_sti_ed=min(ceil((-analwin.win_sti_before+blwin.win_sti_ed)*caframerate),winframe.sti_st-1);
        winframe.blwin_loco_st=max(ceil((-analwin.win_loco_before+blwin.win_loco_st)*caframerate),1);
        winframe.blwin_loco_ed=min(ceil((-analwin.win_loco_before+blwin.win_loco_ed)*caframerate),winframe.loco_st-1);
        
        [soundon,soundoff]=squaredetect(sound_all,0.05);
        [vison,visoff]=squaredetect(diode_all,0.05);
        [airon,airoff]=squaredetect(air_all,0.05);
        [starton,startoff]=squaredetect(startsig_all,0.05);
        
        [caon,caoff]=squaredetect(caframe_all,0.05);
        caon=caon(1:length(caoff));
        
        [pupon,pupoff]=squaredetect(pupilframe_all,0.05);
        pupon=pupon(1:length(pupoff));
        
        [whiskon,whiskoff]=squaredetect(whiskchange_all,0.05);
        whiskon=whiskon(1:length(whiskoff));
        
        [pupilon,pupiloff]=squaredetect(pupilchange_all,0.05);
        pupilon=pupilon(1:length(pupiloff));
        
        [runon,runoff]=squaredetect(wheel_all,0.05);
        runon=runon(1:length(runoff));
        
        [spblon,spbloff]=squaredetect(spontblink_all,0.05);
        spblon=spblon(1:length(spbloff));
        
        class=blinksum_all(:,1);
        [runeachtrial,runspeedeachtrial,pupileachtrial,whiskeachtrial]...
            =deal(NaN(size(class)));
        nTrials=length(class);
        
        % change to time cleanedvector
        spike2tvec=(0:length(air_all))/spike2rate;
        catvec=(caon+caoff)/2/spike2rate;
        puptvec=(pupon+pupoff)/2/spike2rate;
        
        pupframerate=round(max(puptvec)./length(puptvec)*1000);
        % set up event cleanedvectors
        %align to running onset/offset
        runningon=false(length(wheel_all),1);
        i=2;
        while i<=length(runon)
            if (runon(i)+analwin.win_loco_after*spike2rate<length(runningon) ...
                    && runoff(i)+analwin.win_loco_before*spike2rate>0 ...
                    && runoff(i)-runon(i)>analwin.win_loco_after*spike2rate ...
                    && runon(i)-runoff(i-1)>-analwin.win_loco_before*spike2rate)
                runningon((runon(i)+analwin.win_loco_before*spike2rate):(runon(i)+analwin.win_loco_after*spike2rate))=1;
                i=find(runon>runon(i)+analwin.win_loco_after*spike2rate-analwin.win_loco_before*spike2rate,1);
            else
                i=i+1;
            end
        end
        
        runningoff=false(length(wheel_all),1);
        i=1;
        while i<=length(runoff)-1
            if (runoff(i)+analwin.win_loco_after*spike2rate<length(runningoff) ...
                    && runoff(i)+analwin.win_loco_before*spike2rate>0 ...
                    && runoff(i)-runon(i)>analwin.win_loco_after*spike2rate ...
                    && runon(i+1)-runoff(i)>-analwin.win_loco_before*spike2rate)
                runningoff((runoff(i)+analwin.win_loco_before*spike2rate):(runoff(i)+analwin.win_loco_after*spike2rate))=1;
                i=find(runoff>runoff(i)+analwin.win_loco_after*spike2rate-analwin.win_loco_before*spike2rate,1);
            else
                i=i+1;
            end
            
        end
        
        % aligned to cs/us
        cs=zeros(length(diode_all),1);
        for i=1:length(starton)
            cs((vison(i)+analwin.win_sti_before*spike2rate):(vison(i)+analwin.win_sti_after*spike2rate))=1;
            runeachtrial(i)=nanmean(wheel_all((vison(i)-2*spike2rate):vison(i))); % %time spent on 2s before visual cue
            runspeedeachtrial(i)=nanmean(wheelspeed_all((vison(i)-2*spike2rate):vison(i))); % running speed on 2s before visual cue
            st=find(puptvec*spike2rate>vison(i)-2*spike2rate,1);
            ed=find(puptvec*spike2rate>vison(i),1);
            pupileachtrial(i)=nanmean(pupsize_all(st:ed));
            whiskeachtrial(i)=nanmean(whisker_all(st:ed));
        end
        
        % align to spontaneous eyeblinks
        spontbl=false(length(spontblink_all),1);
        i=1;
        while i<=length(spblon)
            if spblon(i)+analwin.win_loco_after*spike2rate<length(spontbl)& spblon(i)+analwin.win_loco_before*spike2rate>0
                spontbl((spblon(i)+analwin.win_loco_before*spike2rate):(spblon(i)+analwin.win_loco_after*spike2rate))=1;
                i=find(spblon>spblon(i)+analwin.win_loco_after*spike2rate-analwin.win_loco_before*spike2rate,1); %exclude after-effects of airpuff
            else
                i=i+1;
            end
        end
        
        % align to whisking
        whiskingon=false(length(whiskchange_all),1);
        i=2;
        while i<=length(whiskon)
            if (whiskon(i)+analwin.win_loco_after*spike2rate<length(whiskingon) ...
                    && whiskon(i)+analwin.win_loco_before*spike2rate>0 ...
                    && whiskoff(i)-whiskon(i)>analwin.win_loco_after*spike2rate ...
                    && whiskon(i)-whiskoff(i-1)>-analwin.win_loco_before*spike2rate)
                whiskingon((whiskon(i)+analwin.win_loco_before*spike2rate):(whiskon(i)+analwin.win_loco_after*spike2rate))=1;
                i=find(whiskon>whiskon(i)+analwin.win_loco_after*spike2rate-analwin.win_loco_before*spike2rate,1);
            else
                i=i+1;
            end
        end
        
        whiskingoff=false(length(whiskchange_all),1);
        i=1;
        while i<=length(whiskoff)-1
            if (whiskoff(i)+analwin.win_loco_after*spike2rate<length(whiskingoff)...
                    && whiskoff(i)+analwin.win_loco_before*spike2rate>0 ...
                    && whiskoff(i)-whiskon(i)>analwin.win_loco_after*spike2rate ...
                    && whiskon(i+1)-whiskoff(i)>-analwin.win_loco_before*spike2rate)
                whiskingoff((whiskoff(i)+analwin.win_loco_before*spike2rate):(whiskoff(i)+analwin.win_loco_after*spike2rate))=1;
                i=find(whiskoff>whiskoff(i)+analwin.win_loco_after*spike2rate-analwin.win_loco_before*spike2rate,1);
            else
                i=i+1;
            end
        end
        
        % align to pupil changes
        pupildila=false(length(pupilchange_all),1);
        i=2;
        while i<=length(pupilon)
            if (pupilon(i)+analwin.win_loco_after*spike2rate<length(pupildila)  ...
                    && pupilon(i)+analwin.win_loco_before*spike2rate>0 ...
                    && pupiloff(i)-pupilon(i)>analwin.win_loco_after*spike2rate ...
                    && pupilon(i)-pupiloff(i-1)>-analwin.win_loco_before*spike2rate)
                pupildila((pupilon(i)+analwin.win_loco_before*spike2rate):(pupilon(i)+analwin.win_loco_after*spike2rate))=1;
                i=find(pupilon>pupilon(i)+analwin.win_loco_after*spike2rate-analwin.win_loco_before*spike2rate,1);
            else
                i=i+1;
            end
        end
        
        pupilcons=false(length(pupilchange_all),1);
        i=1;
        while i<=length(pupilon)-1
            if (pupiloff(i)+analwin.win_loco_after*spike2rate<length(pupildila)...
                    && pupiloff(i)+analwin.win_loco_before*spike2rate>0 ...
                    && pupiloff(i)-pupilon(i)>analwin.win_loco_after*spike2rate ...
                    && pupilon(i+1)-pupiloff(i)>-analwin.win_loco_before*spike2rate)
                pupilcons((pupiloff(i)+analwin.win_loco_before*spike2rate):(pupiloff(i)+analwin.win_loco_after*spike2rate))=1;
                i=find(pupiloff>pupiloff(i)+analwin.win_loco_after*spike2rate-analwin.win_loco_before*spike2rate,1);
            else
                i=i+1;
            end
        end
        
        % pick random duration
        randomduration=false(length(wheel_all),1);
        starton2=[starton;length(wheel_all)];
        for i=1:length(starton2)-1
            st=starton2(i); ed=starton2(i+1);
            offset=(rand(1)*(ed-st)+st);
            offset=round(max(min(offset,ed-(analwin.win_sti_after-analwin.win_sti_before)*spike2rate-10),st+10));
            randomduration(offset:offset+(-analwin.win_sti_before+analwin.win_sti_after)*spike2rate)=1;
        end
        
        
        % get frame indexes for each modality
        for module=1:9
            fprintf(strcat('analyzing module',num2str(module),'\n'));
            switch module
                case 1
                    analyzereg=runningon;
                case 2
                    analyzereg=runningoff;
                case 3
                    analyzereg=cs;
                case 4
                    analyzereg=spontbl;
                case 5
                    analyzereg=whiskingon;
                case 6
                    analyzereg=whiskingoff;
                case 7
                    analyzereg=pupildila;
                case 8
                    analyzereg=pupilcons;
                case 9
                    analyzereg=randomduration;
            end
            
            [eventon,eventoff]=squaredetect(analyzereg,0.05);
            nevents=min(length(eventon),length(eventoff));
            drugindex=NaN(nevents,1);
            if (module~=3 && module~=9)
                frameindex=nan(ceil((-analwin.win_loco_before+analwin.win_loco_after)*caframerate),nevents);
            else
                frameindex=nan(ceil((-analwin.win_sti_before+analwin.win_sti_after)*caframerate),nevents);
            end
            tstart=spike2tvec(diff(analyzereg)>0);
            tend=spike2tvec(diff(analyzereg)<0);
            for ievent=1:nevents
                [spike2points,caframes,pupframes]=getandalign_lt(tstart(ievent),tend(ievent),spike2tvec,catvec,puptvec);
                frameindex(1:length(caframes),ievent)=caframes-caframeoffset;
                drugindex(ievent)=mean(drug_all(spike2points));
            end
            frameindex(end,:)=[];
            
            switch module
                case 1
                    runningon_fi=frameindex;
                    runningon_drug=drugindex;
                case 2
                    runningoff_fi=frameindex;
                    runningoff_drug=drugindex;
                case 3
                    cs_fi=frameindex;
                    cs_drug=drugindex;
                case 4
                    spontbl_fi=frameindex;
                    spontbl_drug=drugindex;
                case 5
                    whiskingon_fi=frameindex;
                    whiskingon_drug=drugindex;
                case 6
                    whiskingoff_fi=frameindex;
                    whiskingoff_drug=drugindex;
                case 7
                    pupildila_fi=frameindex;
                    pupildila_drug=drugindex;
                case 8
                    pupilcons_fi=frameindex;
                    pupilcons_drug=drugindex;
                case 9
                    random_fi=frameindex;
                    random_drug=drugindex;
            end
        end
        
        brainrunningon=NaN([npix,size(runningon_fi)],'single');
        brainrunningoff=NaN([npix,size(runningoff_fi)],'single');
        braincs=NaN([npix,size(cs_fi)],'single');
        brainspontbl=NaN([npix,size(spontbl_fi)],'single');
        brainwhiskingon=NaN([npix,size(whiskingon_fi)],'single');
        brainwhiskingoff=NaN([npix,size(whiskingoff_fi)],'single');
        brainpupildila=NaN([npix,size(pupildila_fi)],'single');
        brainpupilcons=NaN([npix,size(pupilcons_fi)],'single');
        brainrandom=NaN([npix,size(random_fi)],'single');
        
        for iday=2:length(caframeind_split)
            dayind=caframeind_day(iday-1);
            frameind=[caframeind_split(iday-1),caframeind_split(iday)];
            fprintf(strcat('loading dff of day',num2str(dayind),'\n'))
            load(fullfile(pathname,char(animalid(rowind)),'meso_cleaned',...
                strcat(char(animalid(rowind)),'_D',num2str(dayind),'_PixxTime_dff.mat')));
            dffnframe=size(PixxTime_dff,2);
            % search for frame ind in this range, copy brain pixels in
            % to the matrix for each modality
            eventind=(runningon_fi(1,:)>frameind(1) & runningon_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(runningon_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            brainrunningon(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(runningon_fi,1),[]);
            
            eventind=(runningoff_fi(1,:)>frameind(1) & runningoff_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(runningoff_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            brainrunningoff(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(runningoff_fi,1),[]);
            
            eventind=(cs_fi(1,:)>frameind(1) & cs_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(cs_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            braincs(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(cs_fi,1),[]);
            
            eventind=(spontbl_fi(1,:)>frameind(1) & spontbl_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(spontbl_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            brainspontbl(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(spontbl_fi,1),[]);
            
            eventind=(whiskingon_fi(1,:)>frameind(1) & whiskingon_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(whiskingon_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            brainwhiskingon(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(whiskingon_fi,1),[]);
            
            eventind=(whiskingoff_fi(1,:)>frameind(1) & whiskingoff_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(whiskingoff_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            brainwhiskingoff(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(whiskingoff_fi,1),[]);
            
            eventind=(pupildila_fi(1,:)>frameind(1) & pupildila_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(pupildila_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            brainpupildila(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(pupildila_fi,1),[]);
            
            eventind=(pupilcons_fi(1,:)>frameind(1) & pupilcons_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(pupilcons_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            brainpupilcons(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(pupilcons_fi,1),[]);
            
            eventind=(random_fi(1,:)>frameind(1) & random_fi(end,:)<=frameind(2));
            tpframe=caframeind_all(reshape(random_fi(:,eventind),[],1));
            tpframe(tpframe>dffnframe)=dffnframe;
            brainrandom(:,:,eventind)=reshape(PixxTime_dff(allregionspix,tpframe),npix,size(random_fi,1),[]);
        end
        
        % running and pupil aligned around visual stim
        [eventon,eventoff]=squaredetect(cs,0.05);
        nevents=min(length(eventon),length(eventoff));
        tstart=spike2tvec(diff(cs)>0);
        tend=spike2tvec(diff(cs)<0);
        [brainvis_cs,brainair_cs,brainwheel_cs,brainwheelspeed_cs,braincaframe_cs,brainpupilframe_cs]...
            =deal(NaN(spike2rate*(-analwin.win_sti_before+analwin.win_sti_after),nevents));
        
        [brainpupil_cs,brainblink_cs,brainwhisk_cs]...
            =deal(NaN(ceil(pupframerate*(-analwin.win_sti_before+analwin.win_sti_after)),nevents));
        
        for ievent=1:nevents
            [spike2points,caframes,pupframes]=getandalign_lt(tstart(ievent),tend(ievent),spike2tvec,catvec,puptvec);
            vis_cs(1:length(spike2points),ievent)=diode_all(spike2points).*log(max(stipara_all(ievent,1),0.1));
            air_cs(1:length(spike2points),ievent)=air_all(spike2points);
            wheel_cs(1:length(spike2points),ievent)=wheel_all(spike2points);
            wheelspeed_cs(1:length(spike2points),ievent)=wheelspeed_all(spike2points);
            pupilframe_cs(1:length(spike2points),ievent)=pupilframe_all(spike2points);
            caframe_cs(1:length(spike2points),ievent)=caframe_all(spike2points);
            pupil_cs(1:length(pupframes),ievent)=pupsize_all(pupframes);
            blink_cs(1:length(pupframes),ievent)=blink_all(pupframes);
            whisk_cs(1:length(pupframes),ievent)=whisker_all(pupframes);
        end
        
        % blink traced aligned to spontaneous blinks
        [eventon,eventoff]=squaredetect(spontbl,0.05);
        nevents=min(length(eventon),length(eventoff));
        tstart=spike2tvec(diff(spontbl)>0);
        tend=spike2tvec(diff(spontbl)<0);
        blink_spblink=NaN(ceil(pupframerate*(-analwin.win_sti_before+analwin.win_sti_after)),nevents);
        for i=1:nevents
            [spike2points,caframes,pupframes]=getandalign_lt(tstart(i),tend(i),spike2tvec,catvec,puptvec);
            blink_spblink(1:length(pupframes),i)=blink_all(pupframes);
        end
        
        %savefilename=fullfile('C:\Data',char(animalid(rowind)),strcat(char(animalid(rowind)),savesuffix));
        
        save([savefilename,'_caaligned.mat'],...
            'analwin','blwin','winframe','caframerate','spike2rate','pupileachtrial','runeachtrial','runspeedeachtrial','whiskeachtrial',...
            'brainrunningon','brainrunningoff','braincs','brainspontbl','brainrandom','brainwhiskingon','pupildila_drug','pupilcons_drug','cs_drug','spontbl_drug',...
            'brainwhiskingoff','brainpupildila','brainpupilcons','runningon_drug','runningoff_drug','whiskingon_drug','whiskingoff_drug','random_drug',...
            'vis_cs','air_cs','wheel_cs','wheelspeed_cs','pupilframe_cs','caframe_cs','pupil_cs','blink_cs','whisk_cs','blink_spblink','-v7.3');
        
        if type==4 | type==8
            pixmat2tif(mean(braincs,3),allregions,'uint16',200,[savefilename,'_cs_all'],0);
            pixmat2tif(mean(brainrunningon,3),allregions,'uint16',200,[savefilename,'_runningon'],0);
            pixmat2tif(mean(brainrunningoff,3),allregions,'uint16',200,[savefilename,'_runningoff'],0);
            pixmat2tif(mean(brainwhiskingon,3),allregions,'uint16',200,[savefilename,'_whiskingon'],0);
            pixmat2tif(mean(brainwhiskingoff,3),allregions,'uint16',200,[savefilename,'_whiskingoff'],0);
            pixmat2tif(mean(brainpupildila,3),allregions,'uint16',200,[savefilename,'_pupildila'],0);
            pixmat2tif(mean(brainpupilcons,3),allregions,'uint16',200,[savefilename,'_pupilcons'],0);
            pixmat2tif(mean(brainspontbl,3),allregions,'uint16',200,[savefilename,'_spontbl'],0);
            pixmat2tif(mean(brainrandom,3),allregions,'uint16',200,[savefilename,'_random'],0);
        end
        pupthres=nanmedian(pupileachtrial);
        runthres=1;
        whiskthres=0.5;
        if type==4
            daybin=discretize(cell2mat(daysession_all(:,2)),3);
            
            ind=(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==1;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_early_correct'],0);
            ind=(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==1;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_early_incorrect'],0);
            ind=(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==2;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_mid_correct'],0);
            ind=(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==2;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_mid_incorrect'],0);
            ind=(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==3;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_late_correct'],0);
            ind=(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==3;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_late_incorrect'],0);
            
            
            ind=blinksum_all(:,1)==3& stipara_all(:,1)==100 & daybin==1;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_early_catch1'],0);
            ind=blinksum_all(:,1)==2& stipara_all(:,1)==0 & daybin==1;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_early_catch2'],0);
            ind=blinksum_all(:,1)==3& stipara_all(:,1)==100 & daybin==2;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_mid_catch1'],0);
            ind=blinksum_all(:,1)==2& stipara_all(:,1)==0 & daybin==2;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_mid_catch2'],0);
            ind=blinksum_all(:,1)==3& stipara_all(:,1)==100 & daybin==3;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_late_catch1'],0);
            ind=blinksum_all(:,1)==2& stipara_all(:,1)==0 & daybin==3;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_late_catch2'],0);
            
        elseif type==8
            ind=(blinksum_all(:,1)==1|blinksum_all(:,1)==3);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_correct'],0);
            ind=(blinksum_all(:,1)==2|blinksum_all(:,1)==4);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_incorrect'],0);
            ind=(stipara_all(:,1)<10);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_lowcontrast'],0);
            ind=(stipara_all(:,1)>=10 & stipara_all(:,1)<=20);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast'],0);
            ind=(stipara_all(:,1)>20);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_highcontrast'],0);
            
            ind=(stipara_all(:,1)==20 & (blinksum_all(:,1)==1|blinksum_all(:,1)==3));
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast_correct'],0);
            ind=(stipara_all(:,1)==20 & (blinksum_all(:,1)==2|blinksum_all(:,1)==4));
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast_incorrect'],0);
            
            ind=(stipara_all(:,1)==20 & pupileachtrial>=pupthres);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast_largepup'],0);
            ind=(stipara_all(:,1)==20 &  pupileachtrial<pupthres);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast_smallpup'],0);
            ind=(stipara_all(:,1)==20 & whiskeachtrial>=whiskthres);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast_whisking'],0);
            ind=(stipara_all(:,1)==20 &  whiskeachtrial<whiskthres);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast_nowhisking'],0);
            ind=(stipara_all(:,1)==20 & runspeedeachtrial>=runthres);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast_highspeed'],0);
            ind=(stipara_all(:,1)==20 &  runspeedeachtrial<runthres);
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_midcontrast_lowspeed'],0);
            
            daybin=discretize(cell2mat(daysession_all(:,2)),3);
            ind=stipara_all(:,1)==20 & (blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==1;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_early_mid_correct'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==1;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_early_mid_incorrect'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==2;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_mid_mid_correct'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==2;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_mid_mid_incorrect'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==3;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_late_mid_correct'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==3;
            pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_late_mid_incorrect'],0);
            
        elseif type==16
            contrastlvl=unique(stipara_all(:,1));
            for druglvl=1:4
                switch druglvl
                    case 1
                        tit2='_none';
                    case 2
                        tit2='_low';
                    case 3
                        tit2='_mid';
                    case 4
                        tit2='_high';
                end
                
                ind=cs_drug==druglvl;
                pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_all',tit2],0);
                ind=runningon_drug==druglvl;
                pixmat2tif(mean(brainrunningon(:,:,ind),3),allregions,'uint16',200,[savefilename,'_runningon',tit2],0);
                ind=runningoff_drug==druglvl;
                pixmat2tif(mean(brainrunningoff(:,:,ind),3),allregions,'uint16',200,[savefilename,'_runningoff',tit2],0);
                ind=whiskingon_drug==druglvl;
                pixmat2tif(mean(brainwhiskingon(:,:,ind),3),allregions,'uint16',200,[savefilename,'_whiskingon',tit2],0);
                ind=whiskingoff_drug==druglvl;
                pixmat2tif(mean(brainwhiskingoff(:,:,ind),3),allregions,'uint16',200,[savefilename,'_whiskingoff',tit2],0);
                ind=pupildila_drug==druglvl;
                pixmat2tif(mean(brainpupildila(:,:,ind),3),allregions,'uint16',200,[savefilename,'_pupildila',tit2],0);
                ind=pupilcons_drug==druglvl;
                pixmat2tif(mean(brainpupilcons(:,:,ind),3),allregions,'uint16',200,[savefilename,'_pupilcons',tit2],0);
                ind=spontbl_drug==druglvl;
                pixmat2tif(mean(brainspontbl(:,:,ind),3),allregions,'uint16',200,[savefilename,'_spontbl',tit2],0);
                ind=random_drug==druglvl;
                pixmat2tif(mean(brainrandom(:,:,ind),3),allregions,'uint16',200,[savefilename,'_random',tit2],0);
                
                ind=stipara_all(:,17)==druglvl & stipara_all(:,1)==40 & (blinksum_all(:,1)==1|blinksum_all(:,1)==3);
                pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_40_correct',tit2],0);
                ind=stipara_all(:,17)==druglvl & stipara_all(:,1)==40 & (blinksum_all(:,1)==2|blinksum_all(:,1)==4);
                pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_40_incorrect',tit2],0);
                
                for icontrastlvl=1:length(contrastlvl)
                    ind=stipara_all(:,17)==druglvl & stipara_all(:,1)==contrastlvl(icontrastlvl);
                    pixmat2tif(mean(braincs(:,:,ind),3),allregions,'uint16',200,[savefilename,'_cs_',num2str(contrastlvl(icontrastlvl)),tit2],0);
                end
            end
        end
        
        % baseline over x sec window before visual response
        fprintf('Baselining by according window');
        tp=nanmean(brainrunningon(:,winframe.blwin_loco_st:winframe.blwin_loco_ed,:),2);
        brainrunningon_WB=brainrunningon-repmat(tp,1,size(brainrunningon,2),1);
        clear brainrunningon
        tp=nanmean(brainrunningoff(:,winframe.blwin_loco_st:winframe.blwin_loco_ed,:),2);
        brainrunningoff_WB=brainrunningoff-repmat(tp,1,size(brainrunningoff,2),1);
        clear brainrunningoff
        tp=nanmean(braincs(:,winframe.blwin_sti_st:winframe.blwin_sti_ed,:),2);
        braincs_WB=braincs-repmat(tp,1,size(braincs,2),1);
        clear braincs
        tp=nanmean(brainspontbl(:,winframe.blwin_loco_st:winframe.blwin_loco_ed,:),2);
        brainspontbl_WB=brainspontbl-repmat(tp,1,size(brainspontbl,2),1);
        clear brainspontbl
        tp=nanmean(brainwhiskingon(:,winframe.blwin_loco_st:winframe.blwin_loco_ed,:),2);
        brainwhiskingon_WB=brainwhiskingon-repmat(tp,1,size(brainwhiskingon,2),1);
        clear brainwhiskingon
        tp=nanmean(brainwhiskingoff(:,winframe.blwin_loco_st:winframe.blwin_loco_ed,:),2);
        brainwhiskingoff_WB=brainwhiskingoff-repmat(tp,1,size(brainwhiskingoff,2),1);
        clear brainwhiskingoff
        tp=nanmean(brainpupildila(:,winframe.blwin_loco_st:winframe.blwin_loco_ed,:),2);
        brainpupildila_WB=brainpupildila-repmat(tp,1,size(brainpupildila,2),1);
        clear brainpupildila
        tp=nanmean(brainpupilcons(:,winframe.blwin_loco_st:winframe.blwin_loco_ed,:),2);
        brainpupilcons_WB=brainpupilcons-repmat(tp,1,size(brainpupilcons,2),1);
        clear brainpupilcons
        tp=nanmean(brainrandom(:,winframe.blwin_sti_st:winframe.blwin_sti_ed,:),2);
        brainrandom_WB=brainrandom-repmat(tp,1,size(brainrandom,2),1);
        clear brainrandom
        
        save([savefilename,'_caaligned_WB.mat'],...
            'brainrunningon_WB','brainrunningoff_WB','braincs_WB',...
            'brainspontbl_WB','brainrandom_WB','brainwhiskingon_WB',...
            'brainwhiskingoff_WB','brainpupildila_WB','brainpupilcons_WB','-v7.3');
        
        if type==4 | type==8
            pixmat2tif(mean(braincs_WB,3),allregions,'uint16',200,[savefilename,'_WB_cs_all'],0);
            pixmat2tif(mean(brainrunningon_WB,3),allregions,'uint16',200,[savefilename,'_WB_runningon'],0);
            pixmat2tif(mean(brainrunningoff_WB,3),allregions,'uint16',200,[savefilename,'_WB_runningoff'],0);
            pixmat2tif(mean(brainwhiskingon_WB,3),allregions,'uint16',200,[savefilename,'_WB_whiskingon'],0);
            pixmat2tif(mean(brainwhiskingoff_WB,3),allregions,'uint16',200,[savefilename,'_WB_whiskingoff'],0);
            pixmat2tif(mean(brainpupildila_WB,3),allregions,'uint16',200,[savefilename,'_WB_pupildila'],0);
            pixmat2tif(mean(brainpupilcons_WB,3),allregions,'uint16',200,[savefilename,'_WB_pupilcons'],0);
            pixmat2tif(mean(brainspontbl_WB,3),allregions,'uint16',200,[savefilename,'_WB_spontbl'],0);
            pixmat2tif(mean(brainrandom_WB,3),allregions,'uint16',200,[savefilename,'_WB_random'],0);
        end
        if type==4
            daybin=discretize(cell2mat(daysession_all(:,2)),3);
            
            ind=(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==1;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_early_correct'],0);
            ind=(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==1;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_early_incorrect'],0);
            ind=(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==2;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_mid_correct'],0);
            ind=(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==2;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_mid_incorrect'],0);
            ind=(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==3;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_late_correct'],0);
            ind=(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==3;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_late_incorrect'],0);
            
            
            ind=blinksum_all(:,1)==3& stipara_all(:,1)==100 & daybin==1;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_early_catch1'],0);
            ind=blinksum_all(:,1)==2& stipara_all(:,1)==0 & daybin==1;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_early_catch2'],0);
            ind=blinksum_all(:,1)==3& stipara_all(:,1)==100 & daybin==2;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_mid_catch1'],0);
            ind=blinksum_all(:,1)==2& stipara_all(:,1)==0 & daybin==2;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_mid_catch2'],0);
            ind=blinksum_all(:,1)==3& stipara_all(:,1)==100 & daybin==3;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_late_catch1'],0);
            ind=blinksum_all(:,1)==2& stipara_all(:,1)==0 & daybin==3;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_late_catch2'],0);
            
        elseif type==8
            ind=(blinksum_all(:,1)==1|blinksum_all(:,1)==3);
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_correct'],0);
            ind=(blinksum_all(:,1)==2|blinksum_all(:,1)==4);
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_incorrect'],0);
            ind=(stipara_all(:,1)<10);
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_lowcontrast'],0);
            ind=(stipara_all(:,1)>=10 & stipara_all(:,1)<=20);
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_midcontrast'],0);
            ind=(stipara_all(:,1)>20);
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_highcontrast'],0);
            
            ind=(stipara_all(:,1)==20 & (blinksum_all(:,1)==1|blinksum_all(:,1)==3));
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_midcontrast_correct'],0);
            ind=(stipara_all(:,1)==20 & (blinksum_all(:,1)==2|blinksum_all(:,1)==4));
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_midcontrast_incorrect'],0);
            
            daybin=discretize(cell2mat(daysession_all(:,2)),3);
            ind=stipara_all(:,1)==20 & (blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==1;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_early_mid_correct'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==1;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_early_mid_incorrect'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==2;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_mid_mid_correct'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==2;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_mid_mid_incorrect'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==1|blinksum_all(:,1)==3) & daybin==3;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_late_mid_correct'],0);
            ind=stipara_all(:,1)==20 &(blinksum_all(:,1)==2|blinksum_all(:,1)==4) & daybin==3;
            pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_late_mid_incorrect'],0);
        elseif type==16
            contrastlvl=unique(stipara_all(:,1));
            for druglvl=1:4
                switch druglvl
                    case 1
                        tit2='_none';
                    case 2
                        tit2='_low';
                    case 3
                        tit2='_mid';
                    case 4
                        tit2='_high';
                end
                
                ind=cs_drug==druglvl;
                pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_all',tit2],0);
                ind=runningon_drug==druglvl;
                pixmat2tif(mean(brainrunningon_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_runningon',tit2],0);
                ind=runningoff_drug==druglvl;
                pixmat2tif(mean(brainrunningoff_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_runningoff',tit2],0);
                ind=whiskingon_drug==druglvl;
                pixmat2tif(mean(brainwhiskingon_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_whiskingon',tit2],0);
                ind=whiskingoff_drug==druglvl;
                pixmat2tif(mean(brainwhiskingoff_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_whiskingoff',tit2],0);
                ind=pupildila_drug==druglvl;
                pixmat2tif(mean(brainpupildila_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_pupildila',tit2],0);
                ind=pupilcons_drug==druglvl;
                pixmat2tif(mean(brainpupilcons_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_pupilcons',tit2],0);
                ind=spontbl_drug==druglvl;
                pixmat2tif(mean(brainspontbl_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_spontbl',tit2],0);
                ind=random_drug==druglvl;
                pixmat2tif(mean(brainrandom_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_random',tit2],0);
                
                ind=stipara_all(:,17)==druglvl & stipara_all(:,1)==40 & (blinksum_all(:,1)==1|blinksum_all(:,1)==3);
                pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_40_correct',tit2],0);
                ind=stipara_all(:,17)==druglvl & stipara_all(:,1)==40 & (blinksum_all(:,1)==2|blinksum_all(:,1)==4);
                pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_40_incorrect',tit2],0);
                
                for icontrastlvl=1:length(contrastlvl)
                    ind=stipara_all(:,17)==druglvl & stipara_all(:,1)==contrastlvl(icontrastlvl);
                    pixmat2tif(mean(braincs_WB(:,:,ind),3),allregions,'uint16',200,[savefilename,'_WB_cs_',num2str(contrastlvl(icontrastlvl)),tit2],0);
                end
            end
        end
        
        clear brainrunningon_WB brainrunningoff_WB braincs_WB ...
            brainspontbl_WB brainrandom_WB brainwhiskingon_WB ...
            brainwhiskingoff_WB brainpupildila_WB brainpupilcons_WB ...
            vis_cs air_cs wheel_cs wheelspeed_cs pupilframe_cs caframe_cs pupil_cs blink_cs whisk_cs blink_spblink
    end
end

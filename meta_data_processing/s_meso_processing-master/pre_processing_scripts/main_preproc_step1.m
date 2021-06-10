%% output - animal name and suffix (for example, xs_psych_all, etc..)

% Coordinating mesoscopic imaging step1: Concatenate all signals (frame index)
% Crop out marginal regions
% Match timings up with visual cue, pupil and blink traces

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
datapath=strcat('\\128.36.220.173\vivo2\Lan\Meso-imaging');

excel_log = readLogExcel(strcat(datapath,'\animal_log_mesoimaging_lt.xlsx'));
analysisType = 8;

[fieldname, savesuffix, col] = get_preprocessed_files_suffix(analysisType);

%%%%% set up combining criteria here %%%%%
animalnum=find(strcmp(excel_log.Task,'Cond')==1 & strcmp(excel_log.(fieldname),'[]')~=1);
nanimal=length(animalnum);
animalid=excel_log.AnimalID(animalnum,1);

spike2rate=5000;
downsamplefac=20;
if ~isempty(animalnum)
    for rowind=1:nanimal
        caframerate = eval(char(excel_log.Caframerate(animalnum(rowind))));
        if caframerate==10
            caframeimpute=0;
        elseif caframerate>30
            caframeimpute=1;
        end
        adsindex=[]; % records the sequence and combinations of animal/day/session
        
        % generate file list for this animal
        animali=animalnum(rowind);      
        current_animal_name = char(excel_log.AnimalID(animali));
        daytocombine=eval(char(excel_log.(fieldname)(animali)));
        session=1; %%%%%% set session numbers here %%%%%%
        %             temp=combvec(session,daytocombine)';
        temp = [ones(length(daytocombine),1) daytocombine'];
        adsindex=[adsindex;[temp,repmat(rowind,length(session)*length(daytocombine),1)]];
        % col1: session; col2: day (real day); col3: animali;
        filelist=[];
        filelist2=[];
        % generate file list
        for rowindx=1:size(adsindex,1)
            % filename as 'iag_D1'
            %file=fullfile(pathname,current_animal_name,strcat(char(current_animal_name),'_D',num2str(adsindex(rowindx,2)),'_',num2str(adsindex(rowindx,1))));
            file=fullfile(datapath,current_animal_name,strcat(char(current_animal_name),'_D',num2str(adsindex(rowindx,2))));
            fileshort=fullfile(datapath,current_animal_name,'meso_cleaned',strcat(char(current_animal_name),'_D',num2str(adsindex(rowindx,2))));
            filelist=[filelist;file];
            filelist2=[filelist2;fileshort];
        end
        
        % initialize matrices, set to large matrices first
        cacamindex=1;
        caframeind_all=NaN(20000*size(filelist,1),1); %CaRes frame indices over days (20000frames/day *days)
        pupilcamindex=1;
        [blink_all,pupsize_all,whisker_all]=deal(NaN(70000*size(filelist,1),1));
        stipara_all=[];
        blinksum_all=[];
        daysession_all=[];
        
        spike2index=1;
        [pupilframe_all,caframe_all,wheel_all,wheelspeed_all,air_all,...
            diode_all,sound_all,startsig_all,drug_all]=deal(NaN(spike2rate*60*60*size(filelist,1)./downsamplefac,1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load Ca data
        for filenum=1:size(filelist,1) % run through all the files
            filename=char(filelist(filenum,:));
            display (strcat('Type', num2str(analysisType), '_analysis of_', filename));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %crop smr files based on ca camera on and off timing
            fprintf('cropping and concatenating files for %s\n', filename);
            load(strcat(filename,'_binary.mat'));
            load(strcat(filename,'_normpupil.mat')); %norm-pupsize
            load(strcat(filename,'_whisker_clean.mat'));
            load(strcat(filename,'_roi.mat')); %roiint-blink
            load(strcat(filename,'_sti.mat'));
            % fix ITI
            stiparameter(:,15)=[30;stiparameter(1:end-1,15)];
            load(strcat(filename,'_blinksummary.mat'));
            
            %s = whos('-file',strcat(filename2,'_PixxTime_dff.mat'));
            %camax=s.size(2);
            a=dir([filename,'/*.tif']);
            camax=size(a,1);
            %set start and end point
            [starton,startoff]=squaredetect(startsig,0.05);
            excludetrial=find(starton(2:end)-starton(1:end-1)<4*spike2rate)+1; % two trials together, error
            [diodeon,diodeoff]=squaredetect(diode,0.05);
            for itrial=1:length(excludetrial)
                iextrial=excludetrial(itrial);
                startsig(starton(iextrial):startoff(iextrial))=0;
                diode(diodeon(iextrial):diodeoff(iextrial))=0;
                diode(diodeon(iextrial):diodeoff(iextrial))=0;
            end
            stiparameter(excludetrial,:)=[];
            blinksummary(excludetrial,:)=[];
            
            [starton,startoff]=squaredetect(startsig,0.05);
            
            sessionst=starton(1)-10*spike2rate;
            [pupon,pupoff]=squaredetect(pupilframe,0.05);
            pupon(length(roiint)+1:end)=[];pupoff(length(roiint)+1:end)=[];
            
            [caon,caoff]=squaredetect(bluechannel,0.05);
            caon(camax+1:end)=[];caoff(camax+1:end)=[];
            sessioned=min([startoff(end)+10*spike2rate+10,caoff(end)+10,pupoff(end)+10]);
            bluechannel(sessioned(end)+1:end)=0;
            
            usefultrialsstart=find(starton>sessionst,1);
            usefultrialsend=find(starton>sessioned,1);
            if isempty(usefultrialsend)
                usefultrialsend=length(starton);
            else
                startsig(starton(usefultrialsend):end)=0;
                startsig(1:starton(usefultrialsstart)-1)=0;
                air(starton(usefultrialsend):end)=0;
                air(1:starton(usefultrialsstart)-1)=0;
                diode(starton(usefultrialsend):end)=0;
                diode(1:starton(usefultrialsstart)-1)=0;
                sound(starton(usefultrialsend):end)=0;
                sound(1:starton(usefultrialsstart)-1)=0;
                usefultrialsend=usefultrialsend-1;
            end
            
            % concatenate stimulation parameters
            if size(stiparameter,2)<17
                stiparameter(:,17)=0;
            end
            stipara_all=[stipara_all;stiparameter(usefultrialsstart:usefultrialsend,:)];
            blinksum_all=[blinksum_all;blinksummary(usefultrialsstart:usefultrialsend,:)];
            daysession_all=[daysession_all;repmat([animalid(adsindex(filenum,3)),adsindex(filenum,2),adsindex(filenum,1)],usefultrialsend-usefultrialsstart+1,1)];
            
            % coordinate meso frames
            tempdff=find(caon>sessionst & caon<sessioned);
            caframeind_all(cacamindex:cacamindex+size(tempdff,1)-1)=tempdff;
            cacamindex=cacamindex+length(tempdff);
            
            % coordinate pupil frames
            areaii_abs=areaii;
            areaii=areaii_eyenorm;
            
            roiint=sgolayfilt(roiint,1,5);
            baseline=quantile(blinksummary(:,2),0.1);
            maxi=quantile(blinksummary(:,9),0.9);
            roiint(roiint>maxi*1.5)=maxi*1.2;
            roiint=(roiint-baseline)/(maxi-baseline);
            whisk=whisk_binary;
            
            temppup=find(pupon>sessionst & pupon<sessioned)+3;
            pupilframe(1:pupoff(3)+1)=0;
            
            blink_all(pupilcamindex:pupilcamindex+length(temppup)-1)=roiint(temppup);
            pupsize_all(pupilcamindex:pupilcamindex+length(temppup)-1)=areaii(temppup);
            pupsizeabs_all(pupilcamindex:pupilcamindex+length(temppup)-1)=areaii_abs(temppup);
            whisker_all(pupilcamindex:pupilcamindex+length(temppup)-1)=whisk(temppup);
            pupilcamindex=pupilcamindex+length(temppup);
            
            % now crop all else binary data based on calcium camera
            spike2points=length(downsample(sessionst:sessioned,downsamplefac));
            pupilframe_all(spike2index:spike2index+spike2points-1)=downsample(pupilframe(sessionst:sessioned),downsamplefac);
            if caframeimpute==0
                caframe_all(spike2index:spike2index+spike2points-1)=downsample(bluechannel(sessionst:sessioned),downsamplefac);
            elseif caframeimpute==1
                caframe_all(spike2index:spike2index+spike2points-1)=downsample(mesoframe(sessionst:sessioned),downsamplefac);
            end
            wheel_all(spike2index:spike2index+spike2points-1)=downsample(wheel(sessionst:sessioned),downsamplefac);
            wheelspeed_all(spike2index:spike2index+spike2points-1)=downsample(wheelspeed(sessionst:sessioned),downsamplefac);
            air_all(spike2index:spike2index+spike2points-1)=downsample(air(sessionst:sessioned),downsamplefac);
            diode_all(spike2index:spike2index+spike2points-1)=downsample(diode(sessionst:sessioned),downsamplefac);
            sound_all(spike2index:spike2index+spike2points-1)=downsample(sound(sessionst:sessioned),downsamplefac);
            startsig_all(spike2index:spike2index+spike2points-1)=downsample(startsig(sessionst:sessioned),downsamplefac);
            drug_all(spike2index:spike2index+spike2points-1)= stiparameter(1,17)*ones(spike2points,1);
            spike2index=spike2index+spike2points;
        end
        
        pupilframe_all=[0;0;pupilframe_all(3:spike2index-1)];
        caframe_all=[0;0;caframe_all(3:spike2index-1)];
        wheel_all=[0;0;wheel_all(3:spike2index-1)];
        wheelspeed_all=[0;0;wheelspeed_all(3:spike2index-1)];
        air_all=[0;0;air_all(3:spike2index-1)];
        diode_all=[0;0;diode_all(3:spike2index-1)];
        sound_all=[0;0;sound_all(3:spike2index-1)];
        startsig_all=[0;0;startsig_all(3:spike2index-1)];
        
        blink_all=blink_all(1:pupilcamindex-1);
        pupsize_all=pupsize_all(1:pupilcamindex-1);
        whisker_all=whisker_all(1:pupilcamindex-1);
        caframeind_all=caframeind_all(1:cacamindex-1);
        
        % map spontaneous blink timing onto spike2 file
        [pupon,pupoff]=squaredetect(pupilframe_all,0.05);
        [blinkpk_all, blinkloc_all] = findpeaks(blink_all,'MinPeakHeight',0.2,'MinPeakDistance',0.5,'MinPeakProminence',0.2);
        spontblink_all=zeros(size(pupilframe_all));
        
        for i=1:length(blinkloc_all)
            st=max(pupon(blinkloc_all(i))-2*spike2rate/downsamplefac-1,1);
            ed=min(pupon(blinkloc_all(i))+spike2rate/downsamplefac,length(spontblink_all));
            if sum(diode_all(st:ed))==0 % no stimulus during that period
                spontblink_all((st+2*spike2rate/downsamplefac):(st+2*spike2rate/downsamplefac+1))=1;
            end
        end
        
        % map whisking timing onto spike2
        whiskchange_all=zeros(size(diode_all));
        whisker_all(1:2)=0;
        [tpon,tpoff]=squaredetect(whisker_all,0.05);%frame indices
        tpon=tpon(1:length(tpoff));
        for ion=1:length(tpon)
            whiskchange_all(pupon(tpon(ion)):pupoff(tpoff(ion)))=1;
        end
        
        % map pupil dilation timing onto spike2
        pupilchange_all=zeros(size(diode_all));
        [tpon,tpoff]=squaredetect(pupsize_all>0.05,0.05);%frame indices
        if tpon(1)>tpoff(1)
            tpoff(1)=[];
        else
            tpon=tpon(1:length(tpoff));
        end
        for ion=1:min(length(tpon),length(tpoff))
            pupilchange_all(pupon(tpon(ion)):pupoff(tpoff(ion)))=1;
        end
        
        % calculate how much spontaneous blinks per second ITI/ pupilsize/ running/ whisking before each trial
        [spblinkeachtrial,pupileachtrial,pupilabseachtrial,runeachtrial,runspeedeachtrial,whiskeachtrial]...
            =deal(zeros(size(stipara_all,1),1));
        [starton,~]=squaredetect(startsig_all,0.05);
        spontst=round(starton(1)-5*spike2rate/downsamplefac);
        sponted=round(starton(1)-0.1*spike2rate/downsamplefac);
        spblinkeachtrial(1)=sum(spontblink_all(spontst:sponted))/4.9/2;
        for itrial=2:length(starton)
            spontst=starton(itrial-1)+4*spike2rate/downsamplefac;
            sponted=starton(itrial)-0.1*spike2rate/downsamplefac;
            spblinkeachtrial(itrial)=sum(spontblink_all(spontst:sponted))./(stipara_all(itrial,15)-4.1)/2;
        end
        
        [pupon,~]=squaredetect(pupilframe,0.05);
        for itrial=1:length(starton)
            pupst=find(pupon>round(starton(itrial)-2*spike2rate/downsamplefac),1);
            puped=find(pupon>round(starton(itrial)),1)-1;
            pupileachtrial(itrial)=nanmean(pupsize_all(pupst:puped));
            pupilabseachtrial(itrial)=nanmean(pupsizeabs_all(pupst:puped));
            whiskeachtrial(itrial)=nanmean(whisker_all(pupst:puped));
            
            runst=round(starton(itrial)-2*spike2rate/downsamplefac);
            runed=round(starton(itrial))-1;
            runeachtrial(itrial)=nanmean(wheel_all(runst:runed));
            runspeedeachtrial(itrial)=nanmean(wheelspeed_all(runst:runed));
        end
        
        savefilename=char(strcat(datapath,'\',char(current_animal_name),'\',char(current_animal_name),savesuffix,'.mat'));
        fprintf('saving %s\n', savefilename);
        save(savefilename,'spike2rate','caframerate','downsamplefac','blink_all','blinkloc_all','blinkpk_all',...
            'blinksum_all','pupsize_all','pupsizeabs_all','whisker_all','pupilframe_all','spontblink_all','whiskchange_all','pupilchange_all',...
            'caframe_all','wheel_all','wheelspeed_all','air_all','diode_all','sound_all','startsig_all',...
            'stipara_all','daysession_all','caframeind_all','spblinkeachtrial','pupileachtrial','pupilabseachtrial',...
            'runeachtrial','runspeedeachtrial','whiskeachtrial','drug_all','-v7.3');
    end
end

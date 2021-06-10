function detect_timing_spont_blink(animal, daysession, outputpath, overwrite)



rf=0;
for filenum=1:size(daysession,1)
    current_animal=animal(filenum,:);
    data_time_stamp_filename=fullfile(outputpath,current_animal,strcat(current_animal,'_D',num2str(daysession(filenum,1))));
    if exist(strcat(data_time_stamp_filename,'_spontblink.mat'), 'file') && ...
        ~overwrite
        continue;
    end
    display(strcat('spontaneous blink detecting: ',data_time_stamp_filename));
    
    load(strcat(data_time_stamp_filename,'_binary.mat'), 'channels_data');
    load(strcat(data_time_stamp_filename,'_roi.mat'));
    startsig = channels_data.startsig;
    pupilframe = channels_data.pupil_frame;
    [starton,~]=squaredetect(startsig,0.05);
    ntrial=length(starton);
    spontblink=[];
    [pupon,pupoff]=squaredetect(pupilframe,0.05);
    pupon=pupon(1:length(pupoff));
    nframes=length(pupon);
    if rf==0
        [starton,~]=squaredetect(startsig,0.05);
        vison=NaN(ntrial,1);
        visoff=NaN(ntrial,1);
        soundon=NaN(ntrial,1);
        soundoff=NaN(ntrial,1);
        airon=NaN(ntrial,1);
        airoff=NaN(ntrial,1);
        for trialnum=1:ntrial
            [on,off]=squaredetect(diode(starton(trialnum):starton(trialnum)+30000),0.05);
            if ~isempty(on)
                vison(trialnum)=on(1)+starton(trialnum);
                visoff(trialnum)=off(1)+starton(trialnum);
            end
            
            [on,off]=squaredetect(sound(starton(trialnum):starton(trialnum)+30000),0.05);
            if ~isempty(on)
                soundon(trialnum)=on(1)+starton(trialnum);
                soundoff(trialnum)=off(1)+starton(trialnum);
            end
            
            [on,off]=squaredetect(air(starton(trialnum):starton(trialnum)+30000),0.05);
            if ~isempty(on)
                airon(trialnum)=on+starton(trialnum);
                airoff(trialnum)=off+starton(trialnum);
            end
        end
    end
    load(strcat(data_time_stamp_filename,'_roi.mat'), 'roiint');
    roi_sort=sort(roiint);
    blink_bl=mean(roi_sort(1:200));
    blink_pk=mean(roi_sort(end-200:end-100));
    blinkamp=(roiint-blink_bl)./(blink_pk-blink_bl);
    %smoothblinkamp = sgolayfilt(blinkamp,7,21);
    
    cameraoffset=3;
    vtime_pup=timestamp(round((pupon+pupoff)/2));
    vtime_pup=vtime_pup(cameraoffset+1:length(vtime_pup));
    
    % detect spontaneous blinks
    figure('Position',[200 800 1200 200]);
    plot(timestamp,sound,'g');
    hold on;
    plot(timestamp,diode,'r');
    plot(timestamp,air,'k');
    plotlen=min(length(blinkamp),nframes-cameraoffset);
    plot(vtime_pup(1:plotlen),blinkamp(1:plotlen));
    
    if rf==0
        for trial=1:ntrial+1
            spontpk=[];spontloc=[];
            if trial==1
                spontst=10000;
                sponted=starton(trial)-1000;
            elseif trial==ntrial+1
                spontst=starton(trial-1)+20000;
                sponted=pupon(end-5);
            else
                spontst=starton(trial-1)+20000; % 4s after last trial
                sponted=starton(trial)-1000; % 0.2s before vis sti on
            end
            spontpupst=find(vtime_pup>timestamp(spontst),1);
            spontpuped=find(vtime_pup>timestamp(sponted),1);
            
            try
                [spontpk, spontloc] = findpeaks(blinkamp(spontpupst:spontpuped),vtime_pup(spontpupst:spontpuped),'MinPeakHeight',0.1,'MinPeakDistance',0.2,'MinPeakProminence',0.2);
            catch
            end
            spontlocdata=zeros(size(spontpk));
            for t=1:length(spontloc)
                spontlocdata(t)=find(timestamp>spontloc(t),1);
            end
            plot(spontloc,spontpk,'ro');
            % save detection results
            spontblink=[spontblink;[spontloc,spontlocdata,spontpk]]; %#ok<AGROW>
        end
        
    elseif rf==1
        [starton,startoff]=squaredetect(diff(diode)>0,0.05);
        spontst=starton(1);
        sponted=startoff(end);
        spontpupst=find(vtime_pup>timestamp(spontst),1);
        spontpuped=find(vtime_pup>timestamp(sponted),1);
        [spontpk, spontloc] = findpeaks(blinkamp(spontpupst:spontpuped),vtime_pup(spontpupst:spontpuped),'MinPeakHeight',0.1,'MinPeakDistance',0.2,'MinPeakProminence',0.2);
        spontlocdata=zeros(size(spontpk));
        for t=1:length(spontloc)
            spontlocdata(t)=find(timestamp>spontloc(t),1);
        end
        plot(spontloc,spontpk,'ro');
        spontblink=[spontblink;[spontloc,spontlocdata,spontpk]];
    end
    % calculate overall closed time (histogram: time vs eye size)
    %figure;
    %hist(blinkamp(1:plotlen));
    
    save(strcat(data_time_stamp_filename,'_spontblink.mat'),'spontblink');
    
    close all;
    clearvars -except filenum daysession animal pathname rf;
end

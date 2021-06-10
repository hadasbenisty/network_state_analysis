function semi_auto_eye_blink(outputpath, daysession, animal, samplerate, ...
    overwrite)
% type 1 - correct trial
% type 2 - incorrect trial
% type 3 - incorrect response i.e. conditioned blik but no uncondotioned blink (extinction trails, no airpuff)
% type 4 - correct rejection
% type 5 - dumped trial
% circles: green - start of unconditioned response, red - peak of unconditioned response
%          black - start of conditioned response, blue - peak of conditioned response
% if detection is correct click the figure to accept
% if detection is incorrect and trial should be dumped: press any key and then press 1 four times
% if detection is incorrect (either incorrect start and peak locations or incorrect dump trial)
% press any key, this enables mouse cursor and allows 4 clicks.
%       if there is no conditioned response: RIGTH click twice on the figure then left click on start and peak of unconditioned response
%       if there is a conditioned response: LEFT click four times in this order: start(conditioned resp.) - peak(conditioned resp.) - start(unconditioned resp.) - peak(unconditioned resp.)
% after four clicks it will automatically advance to the next trial

for filenum=1:size(daysession,1)
    current_animal=animal(filenum,:);
    data_time_stamp_filename=fullfile(outputpath,current_animal,strcat(current_animal,'_D',num2str(daysession(filenum,1))));
    display(strcat('blink detecting: ',data_time_stamp_filename));
    if exist(strcat(data_time_stamp_filename,'_blink.mat'), 'file') || ...
        ~overwrite
        load(strcat(data_time_stamp_filename,'_blink.mat'),'blink_sti');
    else
        try
        load(strcat(data_time_stamp_filename,'_binary.mat'), 'channels_data');
%         load(strcat(data_time_stamp_filename,'_roi.mat'), 'roiint')
        load(strcat(data_time_stamp_filename,'_sti.mat'), 'stiparameter'); % stimulation
        catch ME
            disp(['file is missing: ' data_time_stamp_filename ]);
            error(ME.message);
        end
        stiparameter=[ones(size(stiparameter(:,1))),stiparameter(:,1:11),zeros(size(stiparameter(:,1))),stiparameter(:,12:size(stiparameter,2))];
        diode = channels_data.diode;
        startsig = zeros(size(diode));
        sound = zeros(size(diode));  % no sound
        %     air = zeros(size(diode));  % no air, use for extinction trials
        startsig(1:end-500)=diode(501:end);
        
        [starton,~]=squaredetect(startsig,0.05);
        ntrial=length(starton);
        
        [pupon,pupoff]=squaredetect(pupilframe,0.05);
        pupon=pupon(1:length(pupoff));
%         nframes=length(pupon);
        
        vison=NaN(ntrial,1);
        visoff=NaN(ntrial,1);
        soundon=NaN(ntrial,1);
        soundoff=NaN(ntrial,1);
        airon=NaN(ntrial,1);
        airoff=NaN(ntrial,1);
        for trialnum=1:ntrial
            [on,off]=squaredetect(diode(starton(trialnum):starton(trialnum)+20000),0.05);
            if ~isempty(on)
                vison(trialnum)=on(1)+starton(trialnum);
                visoff(trialnum)=off(1)+starton(trialnum);
            end
            
            [on,off]=squaredetect(sound(starton(trialnum):starton(trialnum)+20000),0.05);
            if ~isempty(on)
                soundon(trialnum)=on(1)+starton(trialnum);
                soundoff(trialnum)=off(1)+starton(trialnum);
            end
            
            [on,off]=squaredetect(air(starton(trialnum):starton(trialnum)+20000),0.05);
            if ~isempty(on)
                airon(trialnum)=on(1)+starton(trialnum);
                airoff(trialnum)=off(1)+starton(trialnum);
            end
        end
        
        load(strcat(data_time_stamp_filename,'_roi.mat'), 'roiint');
        roi_sort=sort(roiint);
        blink_bl=mean(roi_sort(1:50));
        blink_pk=mean(roi_sort(end-150:end-100));
        blinkamp=(roiint-blink_bl)./(blink_pk-blink_bl);
        %smoothblinkamp = sgolayfilt(blinkamp,7,21);
        
        cameraoffset=3;
        vtime_pup=timestamp(round((pupon+pupoff)/2));
        vtime_pup=vtime_pup(cameraoffset+1:length(vtime_pup));
        
        baseline=NaN(ntrial,1);
        urpeak=NaN(ntrial,1);
        urpeakloc=NaN(ntrial,1);
        uronsetpeak=NaN(ntrial,1);
        uronsetloc=NaN(ntrial,1);
        crpeak=NaN(ntrial,1);
        crpeakloc=NaN(ntrial,1);
        cronsetpeak=NaN(ntrial,1);
        cronsetloc=NaN(ntrial,1);
        class=NaN(ntrial,1);
        blinktrace=NaN(750,ntrial);
        
        % detect peak with in the airpuff window (vis-1000 to vis+2000 ms)
        trial=1;
        while trial<=ntrial
            % find the peak, search backward for cr peak
            try
                if stiparameter(trial,1)==1
                    st=vison(trial)-10000;
                    ed=vison(trial)+15000;
                    mid=vison(trial);
                elseif stiparameter(trial,12)==1
                    st=soundon(trial)-10000;
                    ed=soundon(trial)+15000;
                    mid=soundon(trial);
                end
                
                pupst=find(vtime_pup>timestamp(st),1);
                puped=find(vtime_pup>timestamp(ed),1);
                pupmid=find(vtime_pup>timestamp(mid),1);
                tempblink=blinkamp(pupst:puped);
                blinktrace(1:length(tempblink),trial)=tempblink;
                
                figure(2);
                plot(vtime_pup(pupst:puped),tempblink);
                hold on;
                plot(timestamp(st:ed),sound(st:ed),'g');
                plot(timestamp(st:ed),diode(st:ed),'r');
                plot(timestamp(st:ed),air(st:ed),'k');
                
                % if airpuff exists exists (class 1/2/5)
                if any(air(st:ed))>0
                    %urpeak
                    pupurst=find(vtime_pup>timestamp(airon(trial)),1);
                    pupured=pupurst+5;
                    
                    urpeak(trial)=max(blinkamp(pupurst:pupured));
                    urpeakloc(trial)=vtime_pup(find(blinkamp(pupurst:pupured)==urpeak(trial),1)+pupurst-1);
                    
                    %[urpk, urloc] = findpeaks(blinkamp(pupurst:pupured),vtime_pup(pupurst:pupured),'MinPeakHeight',0.5,'MinPeakProminence',0.02);
                    %urpeak(trial)=urpk(urloc==min(urloc));
                    %urpeakloc(trial)=urloc(urloc==min(urloc));
                    
                    %uronset
%                     pupurmid=find(vtime_pup>urpeakloc(trial),1);
                    
                    uronsetpeak(trial)=min(blinkamp(pupurst-2:pupurst+2));
                    ind=find(blinkamp(pupurst-2:pupurst+2)==uronsetpeak(trial),1);
                    uronsetloc(trial)=vtime_pup(pupurst-2+ind);
                    %[uronstpk, uronstloc] = findpeaks(-blinkamp(pupurst-2:pupurmid-1),vtime_pup(pupurst-2:pupurmid-1),'MinPeakHeight',-1.2);
                    %uronsetpeak(trial)= -uronstpk(uronstloc==min(uronstloc));
                    %uronsetloc(trial)=uronstloc(uronstloc==min(uronstloc));
                    
                    if any(diode(st:ed))>0 || any(sound(st:ed))>0 % if there is cue, thendecide cr
                        if stiparameter(trial,1)==1
                            pupcrst=find(vtime_pup>timestamp(vison(trial)),1);
                        elseif stiparameter(trial,12)==1
                            pupcrst=find(vtime_pup>timestamp(soundon(trial)),1);
                        end
                        pupcred=find(vtime_pup>uronsetloc(trial),1);
                        bl=max(mean(blinkamp(pupcrst-5:pupcrst-2)),mean(blinkamp(pupcrst-15:pupcrst-2)));
                        baseline(trial)=bl;
                        
                        %crpeak
                        crpk=nan;crloc=nan;
                        pk=max(blinkamp(pupcrst:pupcred-2));
                        if pk>bl+0.15*(urpeak(trial)-bl) % if cr peak is large enough
                            crpk=pk;
                            crloc=vtime_pup(find(blinkamp(pupcrst:pupcred-1)==crpk,1)+pupcrst-1);
                        end
                        crpeak(trial)=crpk;
                        crpeakloc(trial)=crloc;
                        %[crpk, crloc] = findpeaks(blinkamp(pupcrst:pupcred-1),vtime_pup(pupcrst:pupcred-1),'MinPeakHeight',0.05+bl,'MinPeakDistance',0.05,'MinPeakProminence',0.02);
                        %crpeak(trial)= crpk(find((crpk==max(crpk))>0,1));
                        %crpeakloc(trial)= crloc(find((crpk==max(crpk))>0,1));
                        
                        %cronset
                        pupcrmid=find(vtime_pup>crpeakloc(trial),1);
                        
                        temp=blinkamp(pupcrst-5:pupcrmid-1);
                        difftemp=diff(temp);
                        ind=find(difftemp>0.05 & temp(1:end-1)>bl,1)-3;
                        if isempty(ind) || isnan(crpeak(trial))
                            class(trial)=2;%no cr reaches threshold
                            crpeak(trial)=nan;
                            crpeakloc(trial)=nan;
                        elseif vtime_pup(pupcrst-5+ind)<timestamp(mid)
                            class(trial)=5; %onset of cr precedes cue
                            urpeak(trial)=nan;
                            urpeakloc(trial)=nan;
                            uronsetpeak(trial)=nan;
                            uronsetloc(trial)=nan;
                        elseif vtime_pup(pupcrst-5+ind)>=timestamp(mid)
                            cronsetloc(trial)=vtime_pup(pupcrst-5+ind);
                            cronsetpeak(trial)=blinkamp(pupcrst-5+ind);
                            class(trial)=1;
                        end
                    else %if there is no cue, only airpuff
                        class(trial)=6;
                    end
                    
                else % if airpuff doesn't exist
                    if any(diode(st:ed))>0 || any(sound(st:ed))>0% decide cr
                        %crpeak
                        if stiparameter(trial,1)==1
                            pupcrst=find(vtime_pup>timestamp(vison(trial)),1);
                        elseif stiparameter(trial,12)==1
                            pupcrst=find(vtime_pup>timestamp(soundon(trial)),1);
                        end
                        
                        pupcred=pupcrst+25;
                        bl=max(mean(blinkamp(pupcrst-5:pupcrst-2)),mean(blinkamp(pupcrst-15:pupcrst-2)));
                        
                        %crpeak
                        pk=max(blinkamp(pupcrst:pupcred-2));
                        if pk>bl+0.15*(1-bl) % if cr peak is large enough
                            crpeak(trial)=pk;
                            crpeakloc(trial)=vtime_pup(find(blinkamp(pupcrst:pupcred-1)==crpeak(trial),1)+pupcrst-1);
                        end
                        
                        %[crpk, crloc] = findpeaks(blinkamp(pupcrst:pupcred-1),vtime_pup(pupcrst:pupcred-1),'MinPeakHeight',0.05+bl,'MinPeakDistance',0.05,'MinPeakProminence',0.02);
                        %crpeak(trial)= crpk(find((crpk==max(crpk))>0,1));
                        %crpeakloc(trial)= crloc(find((crpk==max(crpk))>0,1));
                        
                        %cronset
                        pupcrmid=find(vtime_pup>crpeakloc(trial),1);
                        
                        temp=blinkamp(pupcrst-5:pupcrmid-1);
                        difftemp=diff(temp);
                        ind=find(difftemp>0.05 & temp(1:end-1)>bl,1)-3;
                        
                        if isempty(ind) || isnan(crpeak(trial))
                            class(trial)=4;%no cr reaches threshold
                            crpeak(trial)=nan;
                            crpeakloc(trial)=nan;
                        elseif vtime_pup(pupcrst-5+ind)<timestamp(mid)
                            class(trial)=5; %onset of cr precedes cue
                            crpeak(trial)=nan;
                            crpeakloc(trial)=nan;
                            cronsetpeak(trial)=nan;
                            cronsetloc(trial)=nan;
                        elseif vtime_pup(pupcrst-5+ind)>=timestamp(mid)
                            cronsetloc(trial)=vtime_pup(pupcrst-5+ind);
                            cronsetpeak(trial)=blinkamp(pupcrst-5+ind);
                            class(trial)=3;
                        end
                    else % no cue at all
                        class(trial)=6;
                    end
                end
                
                
                %discard if there is too much baseline variation
                if any(air(st:ed))>0
                    pupur=find(vtime_pup>uronsetloc(trial),1);
                    for move=pupmid-30:pupmid-3
                        moverange=max(blinkamp(move-2:move+3))-blinkamp(min(move-2:move+3));
                        if moverange>(max(blinkamp(pupmid-5:pupur+5))-bl)*0.6 || mean(blinkamp(pupst:pupmid))>0.6*(max(blinkamp(pupmid-5:pupur+5)))
                            class(trial)=5;
                            crpeak(trial)=nan;
                            crpeakloc(trial)=nan;
                            cronsetpeak(trial)=nan;
                            cronsetloc(trial)=nan;
                            urpeak(trial)=nan;
                            urpeakloc(trial)=nan;
                            uronsetpeak(trial)=nan;
                            uronsetloc(trial)=nan;
                        end
                    end
                else
                    for move=pupmid-30:pupmid-3
                        moverange=max(blinkamp(move-2:move+3))-blinkamp(min(move-2:move+3));
                        if (moverange>(max(blinkamp(pupmid-30:pupmid+25))-bl)*0.6 && moverange>0.15) || mean(blinkamp(pupst:pupmid))>0.6
                            class(trial)=5;
                            crpeak(trial)=nan;
                            crpeakloc(trial)=nan;
                            cronsetpeak(trial)=nan;e
                            cronsetloc(trial)=nan;
                            urpeak(trial)=nan;
                            urpeakloc(trial)=nan;
                            uronsetpeak(trial)=nan;
                            uronsetloc(trial)=nan;
                            break;
                        end
                    end
                    
                end
                
                
                if class(trial)==1 || class(trial)==3
                    plot(cronsetloc(trial),cronsetpeak(trial),'ko');
                    plot(crpeakloc(trial),crpeak(trial),'bo');
                end
                if class(trial)==1 ||  class(trial)==2
                    plot(urpeakloc(trial),urpeak(trial),'ro');
                    plot(uronsetloc(trial),uronsetpeak(trial),'go');
                end
                
            catch
                disp('needs manual identification!');
                class(trial)=5;
            end
            hold off;
            display(strcat('Trial=',num2str(trial),'Class=',num2str(class(trial))));
            
            key=waitforbuttonpress; %if not correctly detected, use keyboard key=1
            button=1;
            if key==1
                % now identify blinks
                [x,y,button] = ginput(4); %one cronset, one crpeak, one uronset, one urpeak
                
                cronsetloc(trial)=x(1);
                cronsetpeak(trial)=y(1);
                crpeakloc(trial)=x(2);
                crpeak(trial)=y(2);
                
                uronsetloc(trial)=x(3);
                uronsetpeak(trial)=y(3);
                urpeakloc(trial)=x(4);
                urpeak(trial)=y(4);
                
                
                if ~isnan(airon(trial))
                    if button(1)==1
                        class(trial)=1;
                    elseif button(3)==1
                        class(trial)=2;
                    elseif any(button)==49 % if discard, press 1 for 4 times
                        class(trial)=5; % 5 means discard
                    end
                else % if no airpuff given
                    if button(1)==1 && button(3)==3 %only cr no ur
                        class(trial)=3; % false alarm
                    elseif button(1)==3 && button(3)==3
                        class(trial)=4; % correct rejection
                    elseif any(button)==49 % if discard, press 1 for 4 times
                        class(trial)=5; % 5 means discard
                    end
                end
                
            end
            
            if max(button)>49
                trial=max(trial-2,1);
            end
            trial=trial+1;
        end
        
        blink_sti=NaN(4,2,ntrial);
        for trial=1:ntrial
            if class(trial)==1
                blink_sti(1,1,trial)=cronsetloc(trial);
                blink_sti(1,2,trial)=cronsetpeak(trial);
                blink_sti(2,1,trial)=crpeakloc(trial);
                blink_sti(2,2,trial)=crpeak(trial);
                blink_sti(3,1,trial)=uronsetloc(trial);
                blink_sti(3,2,trial)=uronsetpeak(trial);
                blink_sti(4,1,trial)=urpeakloc(trial);
                blink_sti(4,2,trial)=urpeak(trial);
            elseif class(trial)==2
                blink_sti(3,1,trial)=uronsetloc(trial);
                blink_sti(3,2,trial)=uronsetpeak(trial);
                blink_sti(4,1,trial)=urpeakloc(trial);
                blink_sti(4,2,trial)=urpeak(trial);
            elseif class(trial)==3
            elseif class(trial)==3
                blink_sti(1,1,trial)=cronsetloc(trial);
                blink_sti(1,2,trial)=cronsetpeak(trial);
                blink_sti(2,1,trial)=crpeakloc(trial);
                blink_sti(2,2,trial)=crpeak(trial);
            end
        end
        
        save(strcat(data_time_stamp_filename,'_blink.mat'),'blinkamp','blink_sti','blinktrace');
    end
    ntrial = size(blink_sti, 3);
    
    if ~exist(strcat(data_time_stamp_filename,'_blinksummary.mat'), 'file') || ...
        overwrite
        
        blinksummary=NaN(ntrial,11);
        for trial=1:ntrial
            blinksummary(trial,1)=class(trial);
            if ~isnan(blink_sti(1,1,trial)) % cr exists
                blinksummary(trial,1)=class(trial);
                blinksummary(trial,2)=blink_sti(1,2,trial); %baseline
                if ~isnan(vison(trial))
                    blinksummary(trial,3)=blink_sti(1,1,trial)-vison(trial)./samplerate;%cronset-visonset/soundonset
                    blinksummary(trial,4)=blink_sti(2,1,trial)-vison(trial)./samplerate;%crpeakonset-visonset/soundonset
                elseif ~isnan(soundon(trial))
                    blinksummary(trial,3)=blink_sti(1,1,trial)-soundon(trial)./samplerate;%cronset-visonset/soundonset
                    blinksummary(trial,4)=blink_sti(2,1,trial)-soundon(trial)./samplerate;%crpeakonset-visonset/soundonset
                end
                blinksummary(trial,5)=blink_sti(2,2,trial); %crpeak
                blinksummary(trial,6)=blink_sti(3,2,trial); %crsuspeak
                
            end
            
            if ~isnan(blink_sti(3,1,trial)) % if ur exist
                blinksummary(trial,1)=class(trial);
                blinksummary(trial,2)=baseline(trial); %baseline
                blinksummary(trial,7)=blink_sti(3,1,trial)-airon(trial)./samplerate;%uronset-airpuffonset
                blinksummary(trial,8)=blink_sti(4,1,trial)-airon(trial)./samplerate;%urpeakonset-airpuffonset
                blinksummary(trial,9)=blink_sti(4,2,trial); %urpeak
                blinksummary(trial,10)=(blinksummary(trial,5)-blinksummary(trial,2))./(blinksummary(trial,9)-blinksummary(trial,2)); %crurratio
                blinksummary(trial,11)=(blinksummary(trial,6)-blinksummary(trial,2))./(blinksummary(trial,9)-blinksummary(trial,2));
            else
                blinksummary(trial,2)=baseline(trial); %baseline
                blinksummary(trial,10)=(blinksummary(trial,5)-blinksummary(trial,2))./(1-blinksummary(trial,2)); %cr:full blink ratio
            end
            
        end
        
        correctacceptrate=sum(blinksummary(:,1)==1)./sum(blinksummary(:,1)==1 | blinksummary(:,1)==2);
        correctrejectrate=sum(blinksummary(:,1)==4)./sum(blinksummary(:,1)==3 | blinksummary(:,1)==4);
        trashrate=sum(blinksummary(:,1)==5)./sum(blinksummary(:,1)<6);
        
        display(strcat(data_time_stamp_filename,' correct accept rate=', num2str(correctacceptrate)));
        display(strcat(data_time_stamp_filename,' correct reject rate=', num2str(correctrejectrate)));
        display(strcat(data_time_stamp_filename,' trash rate=', num2str(trashrate)));
        
        save(strcat(data_time_stamp_filename,'_blinksummary.mat'),'blinksummary');
    end
    close all;
    clearvars -except filenum daysession animal pathname;
end
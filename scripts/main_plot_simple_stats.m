function main_plot_simple_stats
animals={'xs','xx','xz','xw','xt','xu'};
%plot_vals_by_state(animals)
%plot_time_spent_by_state_spont(animals)
%plot_learning_sponblink(animals)
%plot_exampletimetraces(animals,1,2,102,110)
end


function plot_vals_by_state(animals)
binsN = 32;
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
for animal_i = 1:length(animals)
    animal=animals{animal_i};
 res(animal_i) = load(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'spon_3states'),...
        'low_pup_q','high_pup_q','high_pup_l',...
        'wheel_low_pup_q', 'pupil_low_pup_q',...
        'face_low_pup_q', 'wheel_high_pup_q', 'pupil_high_pup_q', 'face_high_pup_q', ...
        'wheel_high_pup_l', 'pupil_high_pup_l', 'face_high_pup_l');
end
[parcels_names, parcels_region_labels, finalindex, region_lut] = get_allen_meta_parcels;
CondColors=[0,0,0;0.9290 0.6940 0.1250;1,0,0];
binsN=32;
figure;subplot(3,1,1);
set(gcf,'renderer','Painters')
x=[res.low_pup_q];x=x(finalindex,:);
histogram(x(23,:), binsN,'facecolor',CondColors(1,:),'facealpha',.8,'edgecolor','none'); title('Low Q');
xlim([-100 100]);
x=[res.high_pup_q];x=x(finalindex,:);
subplot(3,1,2);set(gcf,'renderer','Painters');histogram(x(23,:), binsN,'facecolor',CondColors(2,:),'facealpha',.8,'edgecolor','none'); title('High Q');
xlim([-100 100]);
x=[res.high_pup_l];x=x(finalindex,:);
subplot(3,1,3);set(gcf,'renderer','Painters');histogram(x(23,:), binsN,'facecolor',CondColors(3,:),'facealpha',.8,'edgecolor','none'); title('High Loc');
xlim([-100 100]);
xlabel('V1');
suptitle('V1 Distribution');
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\spont_v1_hist_3states');


binsN = linspace(0,60,32);

figure;subplot(3,1,1);set(gcf,'renderer','Painters');
histogram([res.wheel_low_pup_q], binsN,'facecolor',CondColors(1,:),'facealpha',.8,'edgecolor','none');xlim([0 60]); title('Low Q');
subplot(3,1,2);set(gcf,'renderer','Painters');histogram([res.wheel_high_pup_q], binsN,'facecolor',CondColors(2,:),'facealpha',.8,'edgecolor','none'); xlim([0 60]);title('High Q');
subplot(3,1,3);set(gcf,'renderer','Painters');histogram([res.wheel_high_pup_l], binsN,'facecolor',CondColors(3,:),'facealpha',.8,'edgecolor','none');xlim([0 60]);xlabel('Running Speed'); title('High Loc');
xlabel('Running Speed'); suptitle('Running Speed Per State');
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\spont_wheel_hist_3states');

binsN = linspace(0,6e3,32);

figure;subplot(3,1,1);set(gcf,'renderer','Painters');
histogram([res.pupil_low_pup_q], binsN,'facecolor',CondColors(1,:),'facealpha',.8,'edgecolor','none');xlim([0 6e3]);title('Low Q');
subplot(3,1,2);set(gcf,'renderer','Painters');histogram([res.pupil_high_pup_q], binsN,'facecolor',CondColors(2,:),'facealpha',.8,'edgecolor','none');xlim([0 6e3]);title('High Q');
subplot(3,1,3);set(gcf,'renderer','Painters');histogram([res.pupil_high_pup_l], binsN,'facecolor',CondColors(3,:),'facealpha',.8,'edgecolor','none');xlim([0 6e3]);xlabel('Pupil'); title('High Loc');
suptitle('Pupil Per State');
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\spont_pupil_hist_3states');

binsN = linspace(-1e4,1e4,32);

figure;subplot(3,1,1);set(gcf,'renderer','Painters');
histogram([res.face_low_pup_q], binsN,'facecolor',CondColors(1,:),'facealpha',.8,'edgecolor','none');xlim([-1e4 1e4]);title('Low Q');
subplot(3,1,2);set(gcf,'renderer','Painters');histogram([res.face_high_pup_q], binsN,'facecolor',CondColors(2,:),'facealpha',.8,'edgecolor','none');xlim([-1e4 1e4]);title('High Q');
subplot(3,1,3);set(gcf,'renderer','Painters');histogram([res.face_high_pup_l], binsN,'facecolor',CondColors(3,:),'facealpha',.8,'edgecolor','none');xlim([-1e4 1e4]);title('High L');
xlabel('Face Map'); suptitle('Face Map Per State');
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\spont_face_hist_3states');

end
function plot_time_spent_by_state_spont(animals)
for animal_i = 1:length(animals)
    animal=animals{animal_i};
    locomotionperiods=[];quiescenceperiods=[];puplow_on_q=[];
    puphigh_on_q=[];puplow_on_loc=[];puphigh_on_loc=[];facelow_on_q=[];
    facehigh_on_q=[];facelow_on_loc=[];facehigh_on_loc=[];
    [~,days_to_process,lengthdays,midpoint_psych]=animaltodays(animal);
    for di = 1:length(days_to_process)
        disp(di)
        filename = fullfile('X:\Hadas\Meso-imaging\lan\',[animal 'psych'],...
            'spike2Features', [animal num2str(days_to_process(di)) 'arousal_state_ITI.mat']);
        if exist(filename, 'file')
            load(filename, 'running_time_traces');
            locomotionperiods = cat(2, locomotionperiods, running_time_traces.locomotionperiods);
            quiescenceperiods = cat(2, quiescenceperiods, running_time_traces.quiescenceperiods);
            puplow_on_q = cat(2, puplow_on_q, running_time_traces.puplow_on_q);
            puphigh_on_q = cat(2, puphigh_on_q, running_time_traces.puphigh_on_q);
            puplow_on_loc = cat(2, puplow_on_loc, running_time_traces.puplow_on_loc);
            puphigh_on_loc = cat(2, puphigh_on_loc, running_time_traces.puphigh_on_loc);
            if isfield(running_time_traces, 'facelow_on_q')
                facelow_on_q = cat(2, facelow_on_q, running_time_traces.facelow_on_q);
                facehigh_on_q = cat(2, facehigh_on_q, running_time_traces.facehigh_on_q);
                facelow_on_loc = cat(2, facelow_on_loc, running_time_traces.facelow_on_loc);
                facehigh_on_loc = cat(2, facehigh_on_loc, running_time_traces.facehigh_on_loc);
            end
        end
    end
    
    q_and_l = size(locomotionperiods,2)+size(quiescenceperiods,2);
    qp_and_l = size(puplow_on_loc,2)+size(puphigh_on_loc,2)+size(puplow_on_q,2)+size(puphigh_on_q	,2);
    timespend.locomotion(animal_i) = size(locomotionperiods,2)/q_and_l;
    timespend.qu(animal_i) = size(quiescenceperiods,2)/q_and_l;
    timespend.lowpup_locomotion(animal_i) = size(puplow_on_loc,2)/qp_and_l;
    timespend.highpup_locomotion(animal_i) = size(puphigh_on_loc,2)/qp_and_l;
    timespend.lowpup_q(animal_i) = size(puplow_on_q,2)/qp_and_l;
    timespend.highpup_q(animal_i) = size(puphigh_on_q,2)/qp_and_l;
    if exist('facelow_on_loc', 'var')
        qf_and_l = size(facelow_on_loc,2)+size(facehigh_on_loc,2)+size(facelow_on_q,2)+size(facehigh_on_q	,2);
        timespend.lowface_locomotion(animal_i) = size(facelow_on_loc,2)/qf_and_l;
        timespend.highface_locomotion(animal_i) = size(facehigh_on_loc,2)/qf_and_l;
        
        timespend.lowface_q(animal_i) = size(facehigh_on_q,2)/qf_and_l;
        timespend.highface_q(animal_i) = size(facehigh_on_q,2)/qf_and_l;
    else
        timespend.lowface_locomotion(animal_i) = nan;
        timespend.highface_locomotion(animal_i) =  nan;
        timespend.lowface_q(animal_i) = nan;
        timespend.highface_q(animal_i) =  nan;
    end
end
names = fieldnames(timespend);
for ni = 1:length(names)
    M.(names{ni}) = nanmean(timespend.(names{ni}));
    S.(names{ni}) = nanstd(timespend.(names{ni}));
end
figure;
barwitherr([S.locomotion S.qu]/sqrt(length(animals)-1), [M.locomotion M.qu]);
set(gca,'XTickLabel', {'Locomotion','Qu'});
ylabel('Fraction of time');


figure;
CondColors=[0,0,0;0,0,1;0.9290 0.6940 0.1250;1,0,0];
allM=[M.lowpup_q M.lowpup_locomotion M.highpup_q M.highpup_locomotion];
allS=[S.lowpup_q S.lowpup_locomotion S.highpup_q S.highpup_locomotion]./sqrt(length(animals)-1);
subplot(1,1,1)
set(gcf,'renderer','Painters')
hold on
for b = 1:4
    bg=bar(b, allM(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
set(gca,'xtick',1:4)
set(gca,'XTickLabel', {'low pup q','low pup loc' 'high pup q' 'highpup loc'});
h = errorbar(1:4,allM, allS,'LineStyle','none','LineWidth',0.5);title('State');
h.Color='k';
set(h, 'marker', 'none'); 
ylabel('Fraction of time');
%barwitherr([S.lowpup_q S.lowpup_locomotion S.highpup_q S.highpup_locomotion]/sqrt(length(animals)-1),[M.lowpup_q M.lowpup_locomotion M.highpup_q M.highpup_locomotion]);
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\spont_time_spent_4states_bypup');


allM=[M.lowface_q M.lowface_locomotion M.highface_q M.highface_locomotion];
allS=[S.lowface_q S.lowface_locomotion S.highface_q S.highface_locomotion]./sqrt(length(animals)-1);
subplot(1,1,1)
set(gcf,'renderer','Painters')
hold on
for b = 1:4
    bg=bar(b, allM(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
set(gca,'xtick',1:4)
set(gca,'XTickLabel', {'low face q','low face loc' 'high face q' 'high face loc'});
h = errorbar(1:4,allM, allS,'LineStyle','none','LineWidth',0.5);title('State');
h.Color='k';
set(h, 'marker', 'none'); 
ylabel('Fraction of time');
%barwitherr([S.lowpup_q S.lowpup_locomotion S.highpup_q S.highpup_locomotion]/sqrt(length(animals)-1),[M.lowpup_q M.lowpup_locomotion M.highpup_q M.highpup_locomotion]);
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\spont_time_spent_4states_byface');
end

function plot_learning_sponblink(animals)
ispsychometric=false;
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions'));
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude');
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory');
spontblinkrate=nan(6,21);correctacceptrate=nan(6,21);falsealarmrate=nan(6,21);
for animal_id=1:length(animals)
    animal=char(animals(animal_id));
    if ispsychometric
       [~,days_to_process]=animaltodays(animal); 
    else
        days_to_process=animaltodayslearning(animal);
    end
    tic
    for dayy=1:length(unique(days_to_process)) %iterate over psychometric days
        %res = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'), strcat(animal,num2str(days_to_process(dayy)),'imaging_time_traces_global.mat')));
        %trial_labels=res.trialslabels.blinksummary(1:size(res.imaging_time_traces.Allen,3));
        res=load(fullfile(strcat('X:\Lan\Meso-imaging\',animal), strcat(animal,'_D',num2str(days_to_process(dayy)),'_blinksummary.mat')));
        trial_labels=res.blinksummary(:,1);
        spontblink=load(fullfile(strcat('X:\Lan\Meso-imaging\',animal), strcat(animal,'_D',num2str(days_to_process(dayy)),'_spontblink.mat')));
        %spont blink size is the # of spont blinks in the session per that
        %day
        timestamp=load(fullfile(strcat('X:\Lan\Meso-imaging\',animal), strcat(animal,'_D',num2str(days_to_process(dayy)),'_binary.mat')));
        %borrowed from lan
        %timestamp.timestamp is total time in sesion- 4seconds following each trial (cut out
        %from spont blink analysis
        %since time is in seconds, i multiply the whole vector elementwise
        %by 0.45 to get spont blink rate for 450 ms
        [timestamp.starton,~]=squaredetect(timestamp.startsig,0.05);
        spontblinkrate(animal_id,dayy)=size(spontblink.spontblink,1)./(timestamp.timestamp(end-1)-4*length(timestamp.starton)).*0.45; %number of spont blinks scaled with recording duration-exlusion duration
        correctacceptrate(animal_id,dayy)=sum(trial_labels==1)./(sum(trial_labels==1)+sum(trial_labels==2));
        falsealarmrate(animal_id,dayy)=sum(trial_labels==3)./(sum(trial_labels==4)+sum(trial_labels==3));
        %failurerate(animal_id,dayy)=sum(trial_labels==5)./length(sti(:,1));        
    end
    disp(strcat(char(animals(animal_id)),num2str(days_to_process(dayy))));
    toc
    %save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'spon_pupilhigh_pupillow'),'pupilhighalldata','pupilhighalldata_t','pupillowalldata_t','pupillowalldata','days_to_process');
    clearvars -except animal_id animals spontblinkrate correctacceptrate falsealarmrate ispsychometric
end
%error bars and plots averaged across 
%maxdays=size(spontblinkrate,2);
figure;
hold on;
errorbar(1:21,nanmean(spontblinkrate,1),nanstd(spontblinkrate,0,1)./sqrt(length(animals)-1),'r.-');
hold on;
errorbar(1:21,nanmean(correctacceptrate,1),nanstd(correctacceptrate,0,1)./sqrt(length(animals)-1),'k.-');
%hold on;
%errorbar(1:21,nanmean(falsealarmrate,1),nanstd(falsealarmrate,0,1)./sqrt(length(animals)-1),'b.-');
xticks(1:21);
xlim([0 22])
if ispsychometric
xticklabels(-2:19);
xlabel('Day Relative to Psychometric Onset');
title('Population Psychometric Corr Rate vs Spont Blink Rate');
mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\','figure_1','psyc_corr_rate_spontblink'), 'all'); 
else
xticklabels(1:21);
xlabel('Day');
title('Population Learning Corr Rate vs Spont Blink Rate');
mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\','figure_1','learning_corr_rate_spontblink'), 'all');
end
end

function plot_exampletimetraces(animalNames,animal_i,day_i,t1,t2)
%animal_i=1;
% day_i=2;
% t1=102;t2=110;
%animal =xs, dayi=2, t1=102,t2=110
%animal =xs, dayi=2, t1=202,t2=210 also 302 to 210
%animalNames = {'xs'};%'xu'   'xs'    };
addpath(genpath('../parcellation/'));
addpath(genpath('../correlation'));
addpath(genpath('../LSSC-higley-master\LSSC-higley-master'));
addpath(genpath('../../../utils/Questionnaire/'));
addpath(genpath('../pre_processing_scripts'));
addpath(genpath('../pre_processing_scripts'));
addpath('X:\Hadas\Meso-imaging\lan\results\code\Functions')
[~, allen_parcels, maskAllen, maskAllenFrontal, allen_ROI_list] = getParcellsByLansAllansAtlass;
[~,days]=animaltodays(char(animalNames(animal_i)));

fltstr = 'spt';

fsspike2=5e3;
%initialize output and input paths
outfigs = 'X:\Hadas\Meso-imaging\lan';
spike2pth0 = 'X:\Hadas\Meso-imaging\lan\spike2data';
addpath(genpath('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master'));

data_smr_path = 'X:\Lan\Meso-imaging\';
cedpath = 'X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\pre_processing_scripts\utils\CEDS64ML';
switch animalNames{animal_i}
    case {'xu','xv','xt','xs'}
        fsimaing=33;
        delay_filt = 500;
    otherwise
        fsimaing=10;
        delay_filt=150;
end
disp([animalNames{animal_i}  ' '  num2str(days(day_i))]);
datapath = ['X:\Hadas\Meso-imaging\lan\' animalNames{animal_i} 'psych\' fltstr '\'];
spike2pth = fullfile(spike2pth0, animalNames{animal_i});
load(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(days(day_i)) '.mat']),'channels_data',...
    'timing', 't_imaging');
parfile = fullfile(datapath, [animalNames{animal_i} '_' num2str(days(day_i)) '_allen.mat']);
pardataAllan = load(parfile);
if fsimaing < 33
    t_imaging=t_imaging(1:2:end);
end
t_imaging = t_imaging(1:end-delay_filt);

if length(t_imaging) > size(pardataAllan.parcels_time_trace,2)
    t_imaging=t_imaging(1:size(pardataAllan.parcels_time_trace,2));
end
stim_timestamps = timing.stimstart/fsspike2;
stim_timestamps=stim_timestamps(1:75);
parcelsLabels.Allen = allen_parcels.regionNum;
roiLabelsbyAllen.Allen = 1:length(allen_parcels.names);
maskByAllen.Allen = maskAllen;
regionLabel.Allen = allen_parcels.regionNum;
isLeftLabel.Allen = repmat(1:2, 1, 28);
regionLabel.nameslegend = {'Rest','Visual','Parietal','Temp','Aud','R-S','S-S','Motor'};
pupil_data = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
num2str(days(day_i))  '_pupil_clean.mat']);
pupilarea = load(pupil_data,  'areaii');
pupil_Norm=pupilarea.areaii; %get pupil camera timepoints
[timing.pupilstart,timing.pupilend]=squaredetect(channels_data.pupil_frame,.5);
pupil_time=timing.pupilstart/fsspike2;
X_Allen =pardataAllan.parcels_time_trace(:, 1:length(t_imaging));
t_spike2 = 1:length(channels_data.wheel)';
%rmpath(genpath('Z:\Lan\Imaging scripts'))
if size(pupil_Norm,1)==length(timing.pupilstart(1:end))
    timing.pupilstart=timing.pupilstart(1:end);
elseif size(pupil_Norm,1)<length(timing.pupilstart(1:end))
    timing.pupilstart=timing.pupilstart(1:size(pupil_Norm,1));
    timing.pupilend=timing.pupilend(1:size(pupil_Norm,1));
elseif size(pupil_Norm,1)>length(timing.pupilstart(1:end))
    pupil_Norm=pupil_Norm(1:length(timing.pupilstart),:);
end
t_spike2=t_spike2/5000;
firstindx_run=findClosestDouble(t_spike2,t_imaging(1));
firstindx_pup=findClosestDouble(pupil_time,t_imaging(1));
t_spike2(1:firstindx_run)=[];channels_data.wheel(1:firstindx_run)=[];
pupil_time(1:firstindx_pup)=[];pupil_Norm(1:firstindx_pup)=[];
figure;
set(gcf,'renderer','Painters')
ax1=subplot(3,1,1)
plot(t_imaging(t1*fsimaing:t2*fsimaing).',X_Allen(2,t1*fsimaing:t2*fsimaing),'color','k');title('L-V1');
ax2=subplot(3,1,2)
plot(t_spike2(t1*fsspike2:t2*fsspike2),channels_data.wheelspeed(t1*fsspike2:t2*fsspike2).','color','r');title('wheel');
ax3=subplot(3,1,3)
plot(pupil_time(t1*30.3:t2*30.3),pupil_Norm(t1*30.3:t2*30.3).','color',[0.9290 0.6940 0.1250]);title('pupil');
linkaxes([ax1, ax2, ax3],'x');

mysave(gcf,fullfile('X:\Lav\ProcessingDirectory\figure_1',strcat('timetraces_allen',char(animalNames(animal_i)),num2str(days(day_i)),'t',num2str(t1),'to',num2str(t2))),'all');
end
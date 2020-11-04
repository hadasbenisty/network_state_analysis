function supplementary_figs
animals={'xt','xs','xx','xu','xz','xw'};
hemodynamic_comparison(animals)
plotamplitudeplots(animals);
end
function hemodynamic_comparison(animals)
%% for each animal, concatenate imaging data over days, normalize traces, and add trial labels
for cond_i=1:2
    if cond_i==1
        fltstr='spt';
    elseif cond_i==2
        fltstr='pixelwise';
    else
    end
for ir=1:length(animals)
    animal=char(animals(ir));
    [animalDate,dys,~,~]=animaltodays(animal);
    days_to_process=dys;
    CumulativeTrialLabels=[];CumulativeZscored=[];
    for dayy=1:length(unique(days_to_process)) %iterate over psychometric days
        res = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\',fltstr), strcat(animal,num2str(days_to_process(dayy)),'imaging_time_traces_global.mat')));
        %st=findClosestDouble(res.imaging_time_traces.t,-0.5);ed=findClosestDouble(res.imaging_time_traces.t,2);
        zscoredbytrial = zscoretrialwise(res.imaging_time_traces.t, res.imaging_time_traces.Allen, -3, -1); %normalize activity per day
        trials_labels=horzcat(res.trialslabels.injectioncond(1:size(zscoredbytrial,3),:),res.trialslabels.blinksummary(1:size(zscoredbytrial,3),:),res.trialslabels.contrastLabels(1:size(zscoredbytrial,3),:));
        CumulativeTrialLabels = cat(1, CumulativeTrialLabels, trials_labels);   %concatenate trial labels (performance, drug, contrast)
        CumulativeZscored = cat(3, CumulativeZscored, zscoredbytrial); %concatenate normalized trial data       
        clearvars res zscoredbytrial trials_labels
        disp(strcat(animal,num2str(days_to_process(dayy)),fltstr))
    end
    psy2=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\', animalDate,'allen_dataproc_psychometric.mat')));
    AllTrials=CumulativeZscored(:,:,psy2.labels_trials.'<3);
    corr=CumulativeZscored(:,:,psy2.labels_trials.'==1);
    incorr=CumulativeZscored(:,:,psy2.labels_trials.'==2);    
    save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,strcat(fltstr,'trials_hemodynamics')),'AllTrials','corr','incorr');
    clearvars -except ir animals fltstr 
end
end

end
function plotamplitudeplots(animals)
spatialindex=2; %left V1=2 SS-Bl=34
condition='V1';
fields = {'allvector','allvectorcorr','allvectorincorr'}; 
c = cell(length(fields),1);
x = cell2struct(c,fields);
for cond_i=1:2
    if cond_i==1
        fltstr='spt';
    elseif cond_i==2
        fltstr='pixelwise';
    else
    end
for ir=1:length(animals)
    animal=char(animals(ir));
    all_traces=load(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,strcat(fltstr,'trials_hemodynamics')),'AllTrials');
    vector_alltrials=mean(squeeze(mean(all_traces.AllTrials(spatialindex,:,:),1)),2); %correct trials
    vector_corr=mean(squeeze(mean(all_traces.corr(spatialindex,:,:),1)),2); %correct trials
    vector_incorr=mean(squeeze(mean(all_traces.incorr(spatialindex,:,:),1)),2); %correct trials
    t_10=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\','xz','psych\spt'), strcat(animaltodays('xz'),'imaging_time_traces_global.mat')));
    if strcmp(animal,'xt')||strcmp(animal,'xs')||strcmp(animal,'xu')       
        t_33 = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'), strcat(animaltodays(animal),'imaging_time_traces_global.mat')));
        %x_10 = interp1(t_33, x_33, t_10);
        vector_alltrials=interp1(t_33.imaging_time_traces.t, vector_alltrials, t_10.imaging_time_traces.t).';
        vector_corr=interp1(t_33.imaging_time_traces.t, vector_corr, t_10.imaging_time_traces.t).';
        vector_incorr=interp1(t_33.imaging_time_traces.t, vector_incorr, t_10.imaging_time_traces.t).';
    else
    end
    st=findClosestDouble(t_10.imaging_time_traces.t,-2);ed=findClosestDouble(t_10.imaging_time_traces.t,4);
    x.allvector=cat(2,x.allvector,vector_alltrials(st:ed));
    x.allvectorcorr=cat(2,x.allvectorcorr,vector_corr(st:ed));
    x.allvectorincorr=cat(2,x.allvectorincorr,vector_incorr(st:ed));   
    clearvars -except ir animals x t_10 spatialindex condition fltstr cond_i spt pxl
end
    if cond_i==1
        spt=x.allvector;
        sptcorr=x.allvectorcorr;
        sptincorr=x.allvectorincorr;
    elseif cond_i==2
        pxl=x.allvector;
        pxlcorr=x.allvectorcorr;
        pxlincorr=x.allvectorincorr;
    else
    end
end
re=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\','xz','psych\spt'), strcat(animaltodays('xz'),'imaging_time_traces_global.mat')));
stind=findClosestDouble(re.imaging_time_traces.t,-2);enind=findClosestDouble(re.imaging_time_traces.t,4);

figure;
x1 = 0.0; x2 = 0.500;
y1 = -2; y2 = 10;
set(gcf,'renderer','painters');
fill([x1 x1 x2 x2],[y1 y2 y2 y1],[0.8 0.8 0.8],'LineStyle','none')
hold on
x3 = 0.450; x4 = 0.500;
fill([x3 x3 x4 x4],[y1 y2 y2 y1],[0.3020 0.7490 0.9294],'LineStyle','none')
xlim([-2 4])
hold on
shadedErrorBar(transpose(re.imaging_time_traces.t(stind:enind)),mean(spt,2),(std(spt,0,2)./(sqrt(size(spt,2)-1))),'lineprops','g');
hold on
shadedErrorBar(transpose(re.imaging_time_traces.t(stind:enind)),mean(pxl,2),(std(pxl,0,2)./(sqrt(size(pxl,2)-1))),'lineprops','r');
xlabel('Time [sec]');ylabel('Z-DF/F');title(strcat('Avg',{' '},condition,{' '},'Spt vs Pxl'));
hold off
mysave(gcf, fullfile('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude',strcat('Avg','/',condition),'overall_performance_separated'), 'all');


figure;
x1 = 0.0; x2 = 0.500;
y1 = -2; y2 = 10;
set(gcf,'renderer','painters');
fill([x1 x1 x2 x2],[y1 y2 y2 y1],[0.8 0.8 0.8],'LineStyle','none')
hold on
x3 = 0.450; x4 = 0.500;
fill([x3 x3 x4 x4],[y1 y2 y2 y1],[0.3020 0.7490 0.9294],'LineStyle','none')
xlim([-2 4])
hold on
shadedErrorBar(transpose(re.imaging_time_traces.t(stind:enind)),mean(sptcorr,2),(std(sptcorr,0,2)./(sqrt(size(sptcorr,2)-1))),'lineprops','g');
hold on
shadedErrorBar(transpose(re.imaging_time_traces.t(stind:enind)),mean(pxlcorr,2),(std(pxlcorr,0,2)./(sqrt(size(pxlcorr,2)-1))),'lineprops','r');
xlabel('Time [sec]');ylabel('Z-DF/F');title(strcat('Avg',{' '},condition,{' '},'Corr Spt vs Pxl'));
hold off
mysave(gcf, fullfile('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude',strcat('Avg','/',condition),'corr_performance_separated'), 'all');

figure;
x1 = 0.0; x2 = 0.500;
y1 = -2; y2 = 10;
set(gcf,'renderer','painters');
fill([x1 x1 x2 x2],[y1 y2 y2 y1],[0.8 0.8 0.8],'LineStyle','none')
hold on
x3 = 0.450; x4 = 0.500;
fill([x3 x3 x4 x4],[y1 y2 y2 y1],[0.3020 0.7490 0.9294],'LineStyle','none')
xlim([-2 4])
hold on
shadedErrorBar(transpose(re.imaging_time_traces.t(stind:enind)),mean(sptincorr,2),(std(sptincorr,0,2)./(sqrt(size(sptincorr,2)-1))),'lineprops','g');
hold on
shadedErrorBar(transpose(re.imaging_time_traces.t(stind:enind)),mean(pxlincorr,2),(std(pxlincorr,0,2)./(sqrt(size(pxlincorr,2)-1))),'lineprops','r');
xlabel('Time [sec]');ylabel('Z-DF/F');title(strcat('Avg',{' '},condition,{' '},'Incorr Spt vs Pxl'));
hold off
mysave(gcf, fullfile('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude',strcat('Avg','/',condition),'incorr_performance_separated'), 'all');


end
function zscoredbytrial=zscoretrialwise(t,X,st,en)
trialsize=size(X,3);
zscoredbytrial=NaN(56,size(X,2),trialsize);
stind = findClosestDouble(t, st);
enind = findClosestDouble(t, en);
for j=1:trialsize
    zscored=NaN(56,size(X,2)); %initialize vector that will give signal for each parcel (for that trial)
for i=1:56
    trial_alldata=X(i,:,j); %for each trial and each parcel take a signal
    tempmean=nanmean(trial_alldata(stind:enind)); %find mean of baseline 3 seconds before vis stim
    tempstd=nanstd(trial_alldata(stind:enind)); %find sd of baseline 3 seconds before visual stim
    tempdatazscored=(trial_alldata-tempmean)./tempstd; %vectorized equation taking away mean and dividing sd of the trial
    zscored(i,:)=tempdatazscored; 
end
zscoredbytrial(:,:,j)=zscored;
end
end




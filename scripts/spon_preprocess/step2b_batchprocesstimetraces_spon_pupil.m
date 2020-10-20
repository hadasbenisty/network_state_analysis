%%Concatenates spontaneous state “trials” from the previous
%step over all days in psyc testing, reshaped to be parcels over time (for high pupil low pupil).

clear;
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions'));
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude');
animals={'xt','xx','xw','xz','xs','xu'};

cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory');
for ir=1:length(animals)
    animal=char(animals(ir));
    mkdir(strcat('allen_Slope_Amplitude\',animal));
    [~,days_to_process,lengthdays,midpoint_psych]=animaltodays(animal);
    pupilhighalldata=[];pupilhighalldata_t=[];pupillowalldata=[];pupillowalldata_t=[];
    for dayy=1:length(unique(days_to_process)) %iterate over psychometric days
        animalname=strcat(animal,num2str(days_to_process(dayy)));
        if exist(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'arousal_state_ITI.mat')),'file')
            res = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'arousal_state_ITI.mat')));
            %pupil high PT=imaging data, t=time vector
            pupilhighPT=reshape(res.running_time_traces.puphigh_on_q,size(res.running_time_traces.puphigh_on_q,1),[]);
            pupilhight=reshape(res.running_time_traces.t_puphigh_on_q,size(res.running_time_traces.t_puphigh_on_q,1)*size(res.running_time_traces.t_puphigh_on_q,2),1);
            %pupil low
            pupillowPT=reshape(res.running_time_traces.puplow_on_q,size(res.running_time_traces.puplow_on_q,1),[]);
            pupillowt=reshape(res.running_time_traces.t_puplow_on_q,size(res.running_time_traces.t_puplow_on_q,1)*size(res.running_time_traces.t_puplow_on_q,2),1);
            %concatenate across animals
            pupilhighalldata=cat(2,pupilhighalldata,pupilhighPT);
            pupilhighalldata_t=cat(1,pupilhighalldata_t,pupilhight);
            pupillowalldata=cat(2,pupillowalldata,pupillowPT);
            pupillowalldata_t=cat(1,pupillowalldata_t,pupillowt);
        else
        end
        clearvars animalname res
    end
    save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'spon_pupilhigh_pupillow'),'pupilhighalldata','pupilhighalldata_t','pupillowalldata_t','pupillowalldata','days_to_process');
    clearvars -except ir animals
end


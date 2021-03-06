%% Concatenates spontaneous state �trials� from the previous step over all days in psyc testing,
%reshaped to be parcels over time (for running not running).
clear;
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions'));
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude');
animals={'xs','xx','xz','xw','xt','xu'};
%% modified 10/20/2020 to set high pupil in quiescence as not running. the change does not 
% matter much since we're no longer looking at not running as a state
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory');
for ir=1:length(animals)
    animal=char(animals(ir));
    mkdir(strcat('allen_Slope_Amplitude\',animal));
    [~,days_to_process,lengthdays,midpoint_psych]=animaltodays(animal);
    runningalldata=[];runningalldata_t=[];notrunningalldata=[];notrunningalldata_t=[];
    for dayy=1:length(unique(days_to_process)) %iterate over psychometric days
        animalname=strcat(animal,num2str(days_to_process(dayy)));
        if exist(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'arousal_state_ITI.mat')),'file')
            res= load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'arousal_state_ITI.mat')));
            %running data
            runningPT=reshape(res.running_time_traces.puphigh_on_loc,size(res.running_time_traces.puphigh_on_loc,1),[]);
            runningt=reshape(res.running_time_traces.t_puphigh_on_loc,size(res.running_time_traces.t_puphigh_on_loc,1)*size(res.running_time_traces.t_puphigh_on_loc,2),1);
            %not running data
            notrunningPT=reshape(res.running_time_traces.puphigh_on_q,size(res.running_time_traces.puphigh_on_q,1),[]);
            notrunningt=reshape(res.running_time_traces.t_puphigh_on_q,size(res.running_time_traces.t_puphigh_on_q,1)*size(res.running_time_traces.t_puphigh_on_q,2),1);
            %concatenate across animals
            runningalldata=cat(2,runningalldata,runningPT);
            runningalldata_t=cat(1,runningalldata_t,runningt);
            notrunningalldata=cat(2,notrunningalldata,notrunningPT);
            notrunningalldata_t=cat(1,notrunningalldata_t,runningt);
        else
        end
        clearvars animalname res
    end
    save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'spon_running_notrunning'),'runningalldata','runningalldata_t','notrunningalldata_t','notrunningalldata','days_to_process');
    clearvars -except ir animals
end


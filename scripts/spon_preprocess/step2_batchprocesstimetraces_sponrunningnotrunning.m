%% Concatenates spontaneous state “trials” from the previous step over all days in psyc testing,
%reshaped to be parcels over time (for running not running).
clear;
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions'));
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude');
animals={'xs','xx','xz','xw','xt','xu'};

cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory');
for ir=1:length(animals)
    animal=char(animals(ir));
    mkdir(strcat('allen_Slope_Amplitude\',animal));
    [~,days_to_process,lengthdays,midpoint_psych]=animaltodays(animal);
    runningalldata=[];runningalldata_t=[];notrunningalldata=[];notrunningalldata_t=[];
    for dayy=1:length(unique(days_to_process)) %iterate over psychometric days
        animalname=strcat(animal,num2str(days_to_process(dayy)));
        if exist(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'running_ITI.mat')),'file')
            res = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'running_ITI.mat')));
            %running data
            runningPT=reshape(res.running_time_traces.locomotionperiods,size(res.running_time_traces.locomotionperiods,1),[]);
            runningt=reshape(res.running_time_traces.t_wheelon,size(res.running_time_traces.t_wheelon,1)*size(res.running_time_traces.t_wheelon,2),1);
            %not running data
            notrunningPT=reshape(res.running_time_traces.quiescenceperiods,size(res.running_time_traces.quiescenceperiods,1),[]);
            notrunningt=reshape(res.running_time_traces.t_wheeloff,size(res.running_time_traces.t_wheeloff,1)*size(res.running_time_traces.t_wheeloff,2),1);
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


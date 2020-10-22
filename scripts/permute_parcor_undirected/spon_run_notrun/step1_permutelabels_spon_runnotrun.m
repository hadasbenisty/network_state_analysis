clear;
addpath(genpath('X:\Lav\network_state_analysis\functions'));
animals={'xs','xt','xu','xx','xz','xw'};
cd('X:\Lav\network_state_analysis\scripts\permute_parcor_undirected\spon_run_notrun');
%%early vs late timetrace structures
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory');
for ir=1:length(animals)
    animal=char(animals(ir));
    mkdir(strcat('allen_Slope_Amplitude\',animal));
    [~,days_to_process,lengthdays,midpoint_psych]=animaltodays(animal);
    runningalldata=[];runningalldata_t=[];notrunningalldata=[];notrunningalldata_t=[];
    fields = {'locomotionperiods','t_wheelon','quiescenceperiods','t_wheeloff'};
    c = cell(length(fields),1);
    x = cell2struct(c,fields);
    for dayy=1:length(unique(days_to_process)) %iterate over psychometric days
        animalname=strcat(animal,num2str(days_to_process(dayy)));
        if exist(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'running_ITI.mat')),'file')
            res = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'running_ITI.mat')));
            x.locomotionperiods=cat(3,x.locomotionperiods,res.running_time_traces.locomotionperiods);
            x.t_wheelon=cat(2,x.t_wheelon,res.running_time_traces.t_wheelon);
            x.quiescenceperiods=cat(3,x.quiescenceperiods,res.running_time_traces.quiescenceperiods);
            x.t_wheeloff=cat(2,x.t_wheeloff,res.running_time_traces.t_wheeloff);
        else
        end
        clearvars animalname res
    end
    for ki=1:100
        permute_spon_run_lan(x,ki,animal);
    end    
    clearvars -except ir animals
end



clear;
addpath(genpath('X:\Lav\network_state_analysis\functions'));
animals={'xs','xt','xu','xx','xz','xw'};
cd('X:\Lav\network_state_analysis\scripts\permute_parcor_undirected\spon_highpup_lowpup');
for ir=1:length(animals)
    animal=char(animals(ir));
    mkdir(strcat('allen_Slope_Amplitude\',animal));
    [~,days_to_process,lengthdays,midpoint_psych]=animaltodays(animal);
    runningalldata=[];runningalldata_t=[];notrunningalldata=[];notrunningalldata_t=[];
    fields = {'puphigh_on','t_puphigh_on','puplow_on','t_puplow_on'};
    c = cell(length(fields),1);
    x = cell2struct(c,fields);
    for dayy=1:length(unique(days_to_process)) %iterate over psychometric days
        animalname=strcat(animal,num2str(days_to_process(dayy)));
        if exist(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'running_ITI.mat')),'file')
            res = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalname,'running_ITI.mat')));
            x.puphigh_on=cat(3,x.puphigh_on,res.running_time_traces.puphigh_on);
            x.t_puphigh_on=cat(2,x.t_puphigh_on,res.running_time_traces.t_puphigh_on);
            x.puplow_on=cat(3,x.puplow_on,res.running_time_traces.puplow_on);
            x.t_puplow_on=cat(2,x.t_puplow_on,res.running_time_traces.t_puplow_on);
        else
        end
        clearvars animalname res
    end
    for ki=1:100
        permute_spon_pup_lan(x,ki,animal);
    end    
    clearvars -except ir animals
end



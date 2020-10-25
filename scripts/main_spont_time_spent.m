animals={'xx','xz','xs','xw','xu'};

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
barwitherr([S.lowpup_q S.lowpup_locomotion S.highpup_q S.highpup_locomotion]/sqrt(length(animals)-1), ...
    [M.lowpup_q M.lowpup_locomotion M.highpup_q M.highpup_locomotion]);

set(gca,'XTickLabel', {'low pup q','low pup loc' 'high pup q' 'highpup loc'});
ylabel('Fraction of time');
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\spont_time_spent_4states_bypup');

figure;
barwitherr([S.lowface_q S.lowface_locomotion S.highface_q S.highface_locomotion]/sqrt(length(animals)-1), ...
    [M.lowface_q M.lowface_locomotion M.highface_q M.highface_locomotion]);
set(gca,'XTickLabel', {'low face q','low face loc' 'high face q' 'high face loc'});
ylabel('Fraction of time');
mysave(gcf, 'X:\Lav\ProcessingDirectory\parcor_undirected\spont_time_spent_4states_byface');

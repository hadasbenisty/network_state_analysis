function main_plot_simple_stats
animals={'xs','xx','xz','xw','xt','xu'};
plot_vals_by_state(animals)
plot_time_spent_by_state_spont(animals)
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
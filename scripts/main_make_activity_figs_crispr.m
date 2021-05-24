function main_make_activity_figs_crispr

addpath('..\functions');
addpath(genpath('../utils'));
addpath(genpath('../meta_data_processing/'));
procdatapath = 'X:\Hadas\Meso-imaging\CRISPR\analysis_results';

outputfiggolder = 'X:\Hadas\Meso-imaging\CRISPR\Figures\activity';
mkNewDir(outputfiggolder);
statenames_4states = {'low_face','high_face', 'loc',};%,'low_face','high_face'


plot_activity_by_state(outputfiggolder, procdatapath, statenames_4states);
plot_overall_mean_activity_by_state_permutation(outputfiggolder, procdatapath, statenames_4states);
plot_traces_events_averaged_animals(outputfiggolder);
%plot_traces_events(outputfiggolder);
end
function plot_overall_mean_activity_by_state_permutation(outputfiggolder, procdatapath, statenames)

animals_db = get_animals_meta_data_by_csv;
REPS=10;

n=length(animals_db.folder_list);
[~, ~, finalindex] = get_allen_meta_parcels;
dat1=nan(23,n,1);dat2=nan(23,n,1);dat3=nan(23,n,1);
dat1_sh=nan(23,n,REPS);dat2_sha=nan(23,n,REPS);dat2_shb=nan(23,n,REPS);dat3_sh=nan(23,n,REPS);
for ai = 1:n
    if animals_db.toinclude_list(ai) ~= find(strcmp(animals_db.toinclude_lut,'Good'))
        continue;
    end
    animal = animals_db.folder_list{ai};
    if ~isfile(fullfile('X:\Hadas\Meso-imaging\CRISPR\analysis_results\',animal,'arousal_traces_states.mat'))||...
            ~isfile(fullfile('X:\Hadas\Meso-imaging\CRISPR\traces_data', animal, 'Ca_traces_spt_patch11_Allen_dfff.mat'))
        continue;
    end
    res1 = load(fullfile('X:\Hadas\Meso-imaging\CRISPR\analysis_results\',animal,'arousal_traces_states.mat'));
    if ~isfield(res1.segments_arousals, 'high_face')||~isfield(res1.segments_arousals, 'low_face')||...
            isempty(res1.segments_arousals.low_face)||isempty(res1.segments_arousals.high_face)
        continue;
    end
    xa = load(fullfile('X:\Hadas\Meso-imaging\CRISPR\traces_data', animal, 'Ca_traces_spt_patch11_Allen_dfff'));
    parcels_traces=xa.parcels_time_trace(finalindex, :);
    load(fullfile('X:\Hadas\Meso-imaging\CRISPR\traces_data', animal, 'smrx_signals_v4.mat'),'timing');
    t_imaging = timing.bluestart;
    L = min([length(t_imaging) size(parcels_traces,2)]);
    parcels_traces=parcels_traces(:,1:L);
    t_imaging=t_imaging(1:L);
    Ysit = extract_segment(t_imaging, parcels_traces, res1.segments_arousals.sit);
    th = quantile(Ysit',0.005)';
    parcels_traces=parcels_traces-th;
    Yloc = extract_segment(t_imaging, parcels_traces, res1.segments_arousals.loc);
    Yhigh_face = extract_segment(t_imaging, parcels_traces, res1.segments_arousals.high_face);
    Ylow_face = extract_segment(t_imaging, parcels_traces, res1.segments_arousals.low_face);
    labels12 = [zeros(size(Ylow_face,2),1); ones(size(Yhigh_face,2),1)];
    labels23 = [zeros(size(Yhigh_face,2),1); ones(size(Yloc,2),1)];
    Y12 = [Ylow_face  Yhigh_face];
    Y23 = [Yhigh_face  Yloc];
    dat1(:,ai) = nanmean(Y12(:,labels12==0),2);
    dat2(:,ai) = nanmean(Y12(:,labels12==1),2);
    dat3(:,ai) = nanmean(Y23(:,labels23==1),2);
    for r=1:REPS
        
        i=randperm(length(labels12));
        labels12sh=labels12(i);
        i=randperm(length(labels23));
        labels23sh=labels23(i);
        
        dat1_sh(:,ai,r) = nanmean(Y12(:,labels12sh==0),2);
        dat2_sha(:,ai,r) = nanmean(Y12(:,labels12sh==1),2);
        dat2_shb(:,ai,r) = nanmean(Y23(:,labels23sh==0),2);
        dat3_sh(:,ai,r) = nanmean(Y23(:,labels23sh==1),2);
        
    end
    
end




for ti = 1:length(animals_db.type_lut)
    currtype = animals_db.type_lut{ti};
    datM(:,ti,1) = nanmean(dat1(:,animals_db.type_list==ti),2);
    datM(:,ti,2) = nanmean(dat2(:,animals_db.type_list==ti),2);
    datM(:,ti,3) = nanmean(dat3(:,animals_db.type_list==ti),2);
end
L = quantile(datM(:),[0.1, 0.9]);
L=[0.7 2];
figure;l=1;
for si = 1:length(statenames)
for ti = 1:length(animals_db.type_lut)
    
        P = scores_to_heatmap_allen(datM(:, ti,si),0);
        subplot(length(animals_db.type_lut),length(statenames),l);
         plot_vals_heatmap(P,...
                        '',[],  L(1), L(2), 1,colormap(redblue));
                    l=l+1;
                    title([animals_db.type_lut{ti} ' ' statenames{si}]);
    end
    
end
mysave(gcf, fullfile(outputfiggolder, ['mean_activity_by_' num2str(length(statenames)) 'states_per_parcel_heatmap'  ]));

for ti = 1:length(animals_db.type_lut)
    currtype = animals_db.type_lut{ti};
    dat1M(ti) = nanmean(nanmean(dat1(:,animals_db.type_list==ti)));
    dat2M(ti) = nanmean(nanmean(dat2(:,animals_db.type_list==ti)));
    dat3M(ti) = nanmean(nanmean(dat3(:,animals_db.type_list==ti)));
    
    dat1M_sh(ti,:) = nanmean(nanmean(dat1_sh(:,animals_db.type_list==ti,:)));
    dat2M_sha(ti,:) = nanmean(nanmean(dat2_sha(:,animals_db.type_list==ti,:)));
    dat2M_shb(ti,:) = nanmean(nanmean(dat2_shb(:,animals_db.type_list==ti,:)));
    dat3M_sh(ti,:) = nanmean(nanmean(dat3_sh(:,animals_db.type_list==ti,:)));
    
    diff12M(ti) = nanmean(dat2M(:,ti)-dat1M(:,ti));
    diff23M(ti) = nanmean(dat3M(:,ti)-dat2M(:,ti));
    
    diff12M_sh(ti,:) = dat2M_sha(ti,:)-dat1M_sh(ti,:);
    diff23M_sh(ti,:) = dat3M_sh(ti,:)-dat2M_shb(ti,:);
    
    ii = ~isnan((dat1(1,:)))'&animals_db.type_list==ti;
    N(ti)= length(unique(animals_db.animal_list(ii)));
    
    
end
% figure;subplot(2,1,1);bar([dat1M;dat2M;dat3M]');set(gca,'XTickLabel',animals_db.type_lut);legend('low p','high p','loc');
% subplot(2,1,2);bar([mean(dat1M_sh,2) mean(dat2M_sha,2) mean(dat2M_shb,2) mean(dat3M_sh,2)]);set(gca,'XTickLabel',animals_db.type_lut);legend('low p','high p','loc');
figure;
subplot(2,1,1);bar([diff23M' nanmean(diff23M_sh,2)]);set(gca,'XTickLabel',animals_db.type_lut);title('Loc minus High Face');ylabel('\Delta activity');
legend('Data','Shuffled');
subplot(2,1,2);bar([diff12M' nanmean(diff12M_sh,2)]);set(gca,'XTickLabel',animals_db.type_lut);title('High Pupil minus Low Face');ylabel('\Delta activity');
legend('Data','Shuffled');
mysave(gcf,fullfile(outputfiggolder,'overall_activity_low_high_face_loc'));
end


function plot_activity_by_state(outputfiggolder, procdatapath, statenames)
parcels_names = get_allen_meta_parcels;

animals_db = get_animals_meta_data_by_csv;
mean_activity = nan(1, length(statenames), length(animals_db.folder_list));
for i=1:length(animals_db.folder_list)
    
    isimaging_good =animals_db.toinclude_list(i)==3;
    
    if ~isimaging_good
        continue;
    end
    resfile = fullfile(procdatapath,  animals_db.folder_list{i}, 'con_states.mat');
    if isfile(resfile)
        load(resfile);
        for state_i = 1:length(statenames)
            dat = eval(statenames{state_i});
            
            if ~isempty(dat.Allen)
                mean_activity(:, state_i, i) = nanmean(nanmean(dat.Allen, 2));
                mean_activity_parcels(:, state_i, i) = nanmean(dat.Allen, 2);
            end
        end
    end
end
for ti = 1:length(animals_db.type_lut)
    currtype = animals_db.type_lut{ti};
    M(ti,:) = nanmean(mean_activity(:, :, animals_db.type_list==ti), 3);
    S(ti,:) = nanstd(mean_activity(:, :, animals_db.type_list==ti), [],3);
    ii = ~isnan(squeeze(mean_activity(1, :, animals_db.type_list==ti)));
    jj = animals_db.animal_list( animals_db.type_list==ti);
    for kk=1:size(ii,1)
        N(ti,kk)= length(unique(jj(ii(kk,:))));
    end
    M_parcels(:, ti,:) = nanmean(mean_activity_parcels(:, :, animals_db.type_list==ti), 3);
    S_parcels(:, ti,:) = nanstd(mean_activity_parcels(:, :, animals_db.type_list==ti), [],3);
end
for si=1:length(animals_db.type_lut)
    strnamesanimals{si} = [animals_db.type_lut{si} ' N=' num2str(N(si))];
end
for si=1:length(statenames)
    strstatenames{si} = statenames{si};
    strstatenames{si}(strstatenames{si}=='_') = ' ';
end
figure;l=1;
L = quantile(M_parcels(:),[0.1, 0.9]);
for si = 1:length(statenames)
for ti = 1:length(animals_db.type_lut)
    
        P = scores_to_heatmap_allen(M_parcels(:, ti,si),0);
        subplot(length(animals_db.type_lut),length(statenames),l);
         plot_vals_heatmap(P,...
                        '',[],  L(1), L(2), 1,colormap(redblue));
                    l=l+1;
                    title([animals_db.type_lut{ti} ' ' strstatenames{si}]);
    end
    
end
mysave(gcf, fullfile(outputfiggolder, ['mean_activity_by_' num2str(length(statenames)) 'states_per_parcel_heatmap'  ]));

N=min(N');

clrs = get_4states_colors;
if size(M,2) == 2
    clrs = clrs([1 3], :);
end

figure;
for ti = 1:length(animals_db.type_lut)
    subplot(length(animals_db.type_lut),1,ti);
    b=barwitherr(bsxfun(@rdivide, squeeze(S_parcels(:, ti,:)),sqrt(N(ti)-1)), squeeze(M_parcels(:, ti,:)));
    % for ib=1:length(b)
    %     b(ib).FaceColor = clrs(ib,:);
    % end
    set(gca,'XTick',1:length(parcels_names))
    set(gca,'XTickLabel',parcels_names)
    title(strnamesanimals{ti});
    % ylim([-.1 0.15]);
end
mysave(gcf, fullfile(outputfiggolder, ['mean_activity_by_' num2str(length(statenames)) 'states_per_parcel1'  ]));
Mmax=nanmax(M_parcels(:));
Mmin=nanmin(M_parcels(:));
figure;l=1;
for ti = 1:length(animals_db.type_lut)
    for si = 1:length(statenames)
        P = scores_to_heatmap_allen(M_parcels(:,ti,si),0);
        subplot(length(statenames)+1,length(animals_db.type_lut),l);plot_vals_heatmap(P);
        title([strnamesanimals{ti} ' ' strstatenames{si}]);
        set(gca,'CLim',[Mmin Mmax]);c=colorbar;
        l=l+1;
    end
    P = scores_to_heatmap_allen(M_parcels(:,ti,end),0)-scores_to_heatmap_allen(M_parcels(:,ti,1),0);
    subplot(length(statenames)+1,length(animals_db.type_lut),l);plot_vals_heatmap(P);
    title([strnamesanimals{ti} ' ' strstatenames{end} ' minus ' strstatenames{1}]);
    set(gca,'CLim',10*abs(Mmin)*[-1 1]);c=colorbar;
    l=l+1;
end
for si=1:size(M_parcels,3)
    for pii=1:size(M_parcels,1)
        x1 = squeeze(mean_activity_parcels(pii,si, animals_db.type_list==1));
        x2 = squeeze(mean_activity_parcels(pii,si, animals_db.type_list==2));
        x3 = squeeze(mean_activity_parcels(pii,si, animals_db.type_list==3));
        [~,p2(pii,si)]=ttest2(x1,x2,'tail','right');
        [~,p3(pii,si)]=ttest2(x1,x3,'tail','right')
    end
end
figure;subplot(2,2,1);imagesc(p2<0.05);subplot(2,2,2);imagesc(p2>0.95)
subplot(2,2,3);imagesc(p3<0.05);subplot(2,2,4);imagesc(p3>0.95)
figure;
for state_i = 1:length(statenames)
    subplot(length(statenames),1,state_i);
    b=barwitherr(bsxfun(@rdivide, squeeze(S_parcels(:, :,state_i)),sqrt(N-1)),...
        squeeze(M_parcels(:, :,state_i)));
    
    set(gca,'XTick',1:length(parcels_names))
    set(gca,'XTickLabel',parcels_names)
    title(strstatenames{state_i});
    legend(strnamesanimals);
    %     ylim([-0.04 0.14]);
end

mysave(gcf, fullfile(outputfiggolder, ['mean_activity_by_' num2str(length(statenames)) 'states_per_parcel2'  ]));

sem = bsxfun(@rdivide, S,sqrt(N'-1));
figure;b=barwitherr(sem, M);
% for ib=1:length(b)
%     b(ib).FaceColor = clrs(ib,:);
% end
legend(strstatenames);

set(gca,'XTickLabel', strnamesanimals)
ylabel('Mean Activity')

mysave(gcf, fullfile(outputfiggolder, ['mean_activity_by_' num2str(length(statenames)) 'states_populations'  ]));

end
function plot_traces_events_averaged_animals(outputfiggolder)

animals_db = get_animals_meta_data_by_csv;

dffpath = 'X:\Hadas\Meso-imaging\CRISPR\traces_data';
[~, ~, finalindex.Allen] = get_allen_meta_parcels;
validsessions = animals_db.toinclude_list==3;
trials_wheelon = nan(23,71,length(validsessions));
trials_wheeloff= nan(23,71,length(validsessions));
trials_airpuff = nan(23,71,length(validsessions));
trials_high_pupil = nan(23,71,length(validsessions));
for i=1:length(validsessions)
    if ~validsessions(i)
        continue;
    end
    str = animals_db.folder_list{i};
    str(str=='/') = '_';
    str(str=='\') = '_';
    spike2pth = fullfile('X:\Hadas\Meso-imaging\CRISPR\traces_data\', animals_db.folder_list{i});
    if ~isfile(fullfile(spike2pth, 'smrx_signals_v4.mat'))
        continue;
    end
    load(fullfile(spike2pth, 'smrx_signals_v4.mat'), 'timing');
    
    datafile_allen = fullfile(dffpath, animals_db.folder_list{i},  'Ca_traces_spt_patch11_Allen_dfff.mat');
    if ~isfile(datafile_allen)
        continue;
    end
    load(datafile_allen, 'parcels_time_trace');
    parcels_time_trace=parcels_time_trace(finalindex.Allen, :);
    t_imaging = timing.bluestart;
    segmentfile = fullfile('X:\Hadas\Meso-imaging\CRISPR\analysis_results',animals_db.folder_list{i},  'arousal_traces_states.mat');
    if  ~exist(segmentfile,'file')
        continue;
    end
    load(segmentfile, 'segments_arousals');
    Y = extract_segment(t_imaging, parcels_time_trace, segments_arousals.sit);
    Xa = bsxfun(@minus, parcels_time_trace, quantile(Y',.0050)');
    
    
    
    
    before_win=2;
    
    after_win=5;fsimaing=10;
    
    [trials_data, trials_inds] = extract_trials_by_onsets(t_imaging, timing.wheelOff, before_win, after_win, ...
        fsimaing, Xa);
    trials_wheeloff(:,:,i) = mean(trials_data,3);
    
    
    [trials_data, trials_inds] = extract_trials_by_onsets(t_imaging, timing.wheelOn, before_win, after_win, ...
        fsimaing, Xa);
    trials_wheelon(:,:,i) = mean(trials_data,3);
    if isfield(timing, 'airpuffstart') && ~isempty(timing.airpuffstart)
        [trials_data, trials_inds] = extract_trials_by_onsets(t_imaging, timing.airpuffstart, before_win, after_win, ...
            fsimaing, Xa);
        trials_airpuff(:,:,i) = mean(trials_data,3);
    end
    
end

for ci=1:length(animals_db.type_lut)
    M1(:,:,ci) = nanmean(trials_wheelon(:, :,animals_db.type_list==ci),3);
    S1(:,:,ci) = nanstd(trials_wheelon(:, :,animals_db.type_list==ci),[],3);
    
    M2(:,:,ci) = nanmean(trials_airpuff(:,  :,animals_db.type_list==ci),3);
    S2(:,:,ci) = nanstd(trials_airpuff(:,  :,animals_db.type_list==ci),[],3);
    
    M3(:,:,ci) = nanmean(trials_high_pupil(:,  :,animals_db.type_list==ci),3);
    S3(:,:,ci) = nanstd(trials_high_pupil(:,  :,animals_db.type_list==ci),[],3);
    
    M4(:,:,ci) = nanmean(trials_wheeloff(:,  :,animals_db.type_list==ci),3);
    S4(:,:,ci) = nanstd(trials_wheeloff(:,  :,animals_db.type_list==ci),[],3);
    
    ii = animals_db.type_list==ci&~squeeze(isnan(trials_wheelon(1,1,:)));
    nbytype(ci) =  length(unique(animals_db.animal_list(ii)));
end
parcels_names = get_allen_meta_parcels;

inds = [1 14 22];
figure;
for k=1:length(inds)
    subplot(length(inds),1,k);
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M4(inds(k),:,1),S4(inds(k),:,1)/sqrt(nbytype(1)-1));
    hold all;
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M4(inds(k),:,2),S4(inds(k),:,2)/sqrt(nbytype(2)-1),'lineprops','r');
    hold all;
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M4(inds(k),:,3),S4(inds(k),:,3)/sqrt(nbytype(3)-1),'lineprops','b');
    title(parcels_names{inds(k)});
end
c=get(gca,'Children');
legend(c(end:-4:1),animals_db.type_lut);
xlabel('Time [sec]');
suptitle('Activity - Wheel Offset');
mysave(gcf,fullfile(outputfiggolder, 'parcels_wheel_offset'));


figure;
for k=1:length(inds)
    subplot(length(inds),1,k);
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M3(inds(k),:,1),S3(inds(k),:,1)/sqrt(nbytype(1)-1));
    hold all;
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M3(inds(k),:,2),S3(inds(k),:,2)/sqrt(nbytype(2)-1),'lineprops','r');
    hold all;
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M3(inds(k),:,3),S3(inds(k),:,3)/sqrt(nbytype(3)-1),'lineprops','b');
    title(parcels_names{inds(k)});
end
c=get(gca,'Children');
legend(c(end:-4:1),animals_db.type_lut);
xlabel('Time [sec]');
suptitle('Activity - Wheel Onset');
mysave(gcf,fullfile(outputfiggolder, 'parcels_pupil_onset'));



inds = [1 14 22];
figure;
for k=1:length(inds)
    subplot(length(inds),1,k);
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M1(inds(k),:,1),S1(inds(k),:,1)/sqrt(nbytype(1)-1));
    hold all;
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M1(inds(k),:,2),S1(inds(k),:,2)/sqrt(nbytype(2)-1),'lineprops','r');
    hold all;
    shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M1(inds(k),:,3),S1(inds(k),:,3)/sqrt(nbytype(3)-1),'lineprops','b');
    title(parcels_names{inds(k)});
end
c=get(gca,'Children');
legend(c(end:-4:1),animals_db.type_lut);
xlabel('Time [sec]');
suptitle('Activity - Wheel Onset');
mysave(gcf,fullfile(outputfiggolder, 'parcels_wheel_onset'));


for ci=1:length(animals_db.type_lut)
    M(:,ci) = nanmean(nanmean(trials_wheelon(:, :,animals_db.type_list==ci),1),3);
    S(:,ci) = nanstd(nanmean(trials_wheelon(:, :,animals_db.type_list==ci),1),[],3);
    
end
figure;shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,1),S(:,1)/sqrt(nbytype(1)-1));
hold all;
shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,2),S(:,2)/sqrt(nbytype(2)-1),'lineprops','r');
hold all;
shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,3),S(:,3)/sqrt(nbytype(3)-1),'lineprops','b');
c=get(gca,'Children');
legend(c(end:-4:1),animals_db.type_lut);
xlabel('Time [sec]');
title('Overall Average Activity - Wheel Onset');
mysave(gcf,fullfile(outputfiggolder, 'overall_wheel_onset'));

for ci=1:length(animals_db.type_lut)
    M(:,ci) = nanmean(nanmean(trials_wheeloff(:, :,animals_db.type_list==ci),1),3);
    S(:,ci) = nanstd(nanmean(trials_wheeloff(:, :,animals_db.type_list==ci),1),[],3);
    
end
figure;shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,1),S(:,1)/sqrt(nbytype(1)-1));
hold all;
shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,2),S(:,2)/sqrt(nbytype(2)-1),'lineprops','r');
hold all;
shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,3),S(:,3)/sqrt(nbytype(3)-1),'lineprops','b');
c=get(gca,'Children');
legend(c(end:-4:1),animals_db.type_lut);
xlabel('Time [sec]');
title('Overall Average Activity - Wheel Offset');
mysave(gcf,fullfile(outputfiggolder, 'overall_wheel_offset'));



for ci=1:length(animals_db.type_lut)
    M(:,ci) = nanmean(nanmean(trials_high_pupil(:, :,animals_db.type_list==ci),1),3);
    S(:,ci) = nanstd(nanmean(trials_high_pupil(:, :,animals_db.type_list==ci),1),[],3);
    
end
figure;shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,1),S(:,1)/sqrt(nbytype(1)-1));
hold all;
shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,2),S(:,2)/sqrt(nbytype(2)-1),'lineprops','r');
hold all;
shadedErrorBar(linspace(-before_win, after_win, size(M,1)), M(:,3),S(:,3)/sqrt(nbytype(3)-1),'lineprops','b');
c=get(gca,'Children');
legend(c(end:-4:1),animals_db.type_lut);
xlabel('Time [sec]');
title('Overall Average Activity - Wheel Onset');
mysave(gcf,fullfile(outputfiggolder, 'overall_wheel_onset'));

%
% figure;
% for k=1:length(inds)
%  subplot(length(inds),1,k);
% shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M2(inds(k),:,1),S2(inds(k),:,1)/sqrt(nbytype(1)-1));
% hold all;
% shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M2(inds(k),:,2),S2(inds(k),:,2)/sqrt(nbytype(2)-1),'lineprops','r');
% hold all;
% shadedErrorBar(linspace(-before_win, after_win, size(M1,2)), M2(inds(k),:,3),S2(inds(k),:,3)/sqrt(nbytype(3)-1),'lineprops','b');
% title(parcels_names{inds(k)});
% end
% c=get(gca,'Children');
% legend(c(end:-4:1),animals_db.type_lut)


end
function plot_traces_events(outputfiggolder)

animals_db = get_animals_meta_data_by_csv;

dffpath = 'X:\Hadas\Meso-imaging\CRISPR\traces_data';
[~, ~, finalindex.Allen] = get_allen_meta_parcels;
validsessions = animals_db.toinclude_list==3;
for i=1:length(validsessions)
    if ~validsessions(i)
        continue;
    end
    str = animals_db.folder_list{i};
    str(str=='/') = '_';
    str(str=='\') = '_';
    spike2pth = fullfile('X:\Hadas\Meso-imaging\CRISPR\traces_data\', animals_db.folder_list{i});
    if ~isfile(fullfile(spike2pth, 'smrx_signals_v4.mat'))
        continue;
    end
    load(fullfile(spike2pth, 'smrx_signals_v4.mat'), 'timing');
    
    datafile_allen = fullfile(dffpath, animals_db.folder_list{i},  'Ca_traces_spt_patch11_Allen_dfff.mat');
    if ~isfile(datafile_allen)
        continue;
    end
    load(datafile_allen, 'parcels_time_trace');
    parcels_time_trace=parcels_time_trace(finalindex.Allen, :);
    t_imaging = timing.bluestart;
    if isfield(timing, 'airpuffstart') && ~isempty(timing.airpuffstart)
        plotbyevent(timing.airpuffstart, parcels_time_trace, t_imaging);
        suptitle('AirPuff');
        mysave(gcf,fullfile(outputfiggolder, [str '_AirPuff']));
        
    end
    plotbyevent(timing.wheelOn, parcels_time_trace, t_imaging);
    suptitle('Wheel on');
    mysave(gcf,fullfile(outputfiggolder, [str '_WheelOn']));
    close all;
end
end
function plotbyevent(eventtimes, parcels_time_trace, t_imaging)
N=length(eventtimes);
win = 50;X=nan(23,101,N);
for k=1:N
    ind =  findClosestDouble(t_imaging, eventtimes(k));
    if ind+win > size(parcels_time_trace,2)
        break;
    end
    X(:,:,k) = parcels_time_trace(:,ind-win:ind+win);
    
end
tt=linspace(-5,5,101);
figure;
shadedErrorBar(tt,nanmean(mean(X(:,:,:)),3), nanstd(mean(X(:,:,:)),[],3)/sqrt(N-1));
xlabel('Time [sec]');
ylabel('\Delta F/F');

end
function main_make_activity_figs_crispr

addpath('..\functions');
addpath(genpath('../utils'));
addpath(genpath('../meta_data_processing/'));
procdatapath = 'X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr';

outputfiggolder = 'X:\Hadas\Mesoimaging\crispr\meso_results\figs_crispr\activity';
mkNewDir(outputfiggolder);


plot_activity_by_state(outputfiggolder, procdatapath, {'qui', 'loc'});
plot_activity_by_state(outputfiggolder, procdatapath, {'low_pup_q', 'high_pup_q', 'high_pup_l'});

plot_traces_events(outputfiggolder);
end


function plot_activity_by_state(outputfiggolder, procdatapath, statenames)
if length(statenames) == 2
    datafile = 'spont_data_2states_dfff.mat';
else
    datafile = 'spont_data_3states_dfff.mat';
end
animals_db = get_animals_meta_data_by_csv;
mean_activity = nan(1, length(statenames), length(animals_db.folder_list));
for i=1:length(animals_db.folder_list)
    ispup_good =ispupilgood(animals_db, i);
    isimaging_good =isimaginggood(animals_db, i)&animals_db.toinclude_list(i)==3;
    if length(statenames) == 2
        isvalidsession = isimaging_good;
    else
        isvalidsession = isimaging_good & ispup_good;
    end
    if ~isvalidsession
        continue;
    end
    resfile = fullfile(procdatapath,  animals_db.folder_list{i}, datafile);
    if isfile(resfile)
        load(resfile);
        for state_i = 1:length(statenames)
            dat = eval(statenames{state_i});
            if ~isempty(dat.Allen)
                mean_activity(:, state_i, i) = nanmean(nanmean(dat.Allen, 2));
            end
        end
    end
end
for ti = 1:length(animals_db.type_lut)
    currtype = animals_db.type_lut{ti};
    M(ti,:) = nanmean(mean_activity(:, :, animals_db.type_list==ti), 3);
    S(ti,:) = nanstd(mean_activity(:, :, animals_db.type_list==ti), [],3);
    N(ti,:) = sum(~isnan(mean_activity(1, :, animals_db.type_list==ti)), 3);
end
clrs = get_3states_colors;
if size(M,2) == 2
    clrs = clrs([1 3], :);
end

figure;b=barwitherr(S./sqrt(N-1), M);
for ib=1:length(b)
    b(ib).FaceColor = clrs(ib,:);
end
legend(statenames);
set(gca,'XTickLabel', animals_db.type_lut)
ylabel('Mean \Delta F/F')

mysave(gcf, fullfile(outputfiggolder, ['mean_activity_by_' num2str(length(statenames)) 'states_populations'  ]));

end
function plot_traces_events(outputfiggolder)

animals_db = get_animals_meta_data_by_csv;

dffpath = 'X:\Hadas\Meso-imaging\Antara\final_preprocess\alldata';
[~, ~, finalindex.Allen] = get_allen_meta_parcels;
validsessions = animals_db.isgoodpupil_list'==find(strcmp(animals_db.isgoodpupil_lut,'GOOD'));
for i=1:length(validsessions)
    if ~validsessions(i)
        continue;
    end
    str = animals_db.folder_list{i};
    str(str=='/') = '_';
    str(str=='\') = '_';
    spike2pth = fullfile('X:\Hadas\Meso-imaging\Antara\data\Antara\AnalyzedData\', animals_db.folder_list{i});
    if ~isfile(fullfile(spike2pth, 'smrx_signals_v3.mat'))
        continue;
    end
    load(fullfile(spike2pth, 'smrx_signals_v3.mat'), 'timing');
    
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